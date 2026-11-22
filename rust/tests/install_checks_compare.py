#!/usr/bin/env python3
"""Compare the installer checks in priva with the inline Python blocks they replaced.

  python install_checks_compare.py check --repo CHECKOUT --binary PRIVA --work-root NEW_DIRECTORY

The reference blocks are CHECKOUT/tests/fixtures/install-checks-8070d3d/*.py,
copied verbatim from scripts/install_utils.sh at commit 8070d3d, where they ran
as `python - ARGS <<'PY'` heredocs:

| Block | Shell function | priva command |
| --- | --- | --- |
| cds_fasta.py | validate_cds_fasta | cds-fasta validate |
| loftee_resources.py | validate_loftee_resources | loftee-resources validate |
| mechanism_evidence.py | validate_gene_pathogenic_mechanism_evidence | mechanism-evidence validate |
| mondo_scope.py | validate_mondo_hpo_scope_cache | mondo-scope validate |
| spliceai2_check.py, spliceai2_register.py | SpliceAI2_install | spliceai2-install check, register |
| vep_release_pins.py | conda_install_vep | vep-release pins |

Each case runs both on its own copy of the same small generated input and
compares: success, stdout, the reason of a failure (the last stderr line,
without the `priva: ` prefix or a traceback's exception class), and for the
commands that write a file (register, pins) the bytes written. Paths are
compared with each copy's folder replaced by ROOT. A case marked `message`
False compares success only: there Python stops with a traceback (KeyError,
AttributeError, UnicodeDecodeError) whose text the Rust command does not copy.
No network, no real resource.
"""
from __future__ import annotations

import argparse
import fcntl
import gzip
import hashlib
import json
import os
from pathlib import Path
import shutil
import sqlite3
import subprocess
import sys

REVISION = "cd1b8b82bcef9c97ac3b3b9a24f0de3a5bbf77f9"


def sha(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def gz(path: Path, text: bytes) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    with gzip.open(path, "wb") as stream:
        stream.write(text)
    return path


def reason(stderr: str, root: str) -> str:
    lines = [line for line in stderr.replace(root, "ROOT").splitlines() if line.strip()]
    if not lines:
        return ""
    last = lines[-1]
    if last.startswith("priva: "):
        return last[len("priva: "):]
    if any(line.startswith("Traceback") for line in lines):
        name, _, message = last.partition(": ")
        if name.isidentifier():
            return message
    return last


class Compare:
    def __init__(self, repo: Path, binary: Path, work: Path):
        self.blocks = repo / "tests/fixtures/install-checks-8070d3d"
        self.binary = binary
        self.work = work
        self.rows = []
        self.number = 0

    def folder(self, side: str) -> Path:
        path = self.work / f"{self.number:03d}" / side
        path.mkdir(parents=True)
        return path

    def run(self, argv, root: Path, env=None):
        result = subprocess.run([str(a) for a in argv], capture_output=True, text=True, env=env)
        return dict(
            success=result.returncode == 0,
            returncode=result.returncode,
            stdout=result.stdout.replace(str(root), "ROOT"),
            reason=reason(result.stderr, str(root)),
            stderr=result.stderr[-2000:],
        )

    def case(self, name, build, reference, native, *, message=True, written=None, hold_lock=None):
        """build(folder) -> context; reference(ctx) and native(ctx) -> argv lists."""
        self.number += 1
        results = {}
        for side, argv_of in (("reference", reference), ("native", native)):
            root = self.folder(side)
            context = build(root)
            holder = None
            if hold_lock:
                holder = open(hold_lock(context), "a")
                fcntl.flock(holder, fcntl.LOCK_EX | fcntl.LOCK_NB)
            try:
                result = self.run(argv_of(context), root)
            finally:
                if holder:
                    holder.close()
            if written:
                path = written(context)
                result["written"] = path.read_bytes().decode("utf-8", "replace").replace(str(root), "ROOT") if path.exists() else None
            results[side] = result
        ref, nat = results["reference"], results["native"]
        passed = ref["success"] == nat["success"]
        if ref["success"]:
            passed = passed and ref["stdout"] == nat["stdout"]
        elif message:
            passed = passed and ref["reason"] == nat["reason"]
        if written:
            passed = passed and ref["written"] == nat["written"]
        self.rows.append(dict(name=name, passed=passed, message_checked=message and not ref["success"], **results))
        status = "ok " if passed else "MISMATCH"
        print(f"{status} {self.number:03d} {name}: {'success' if ref['success'] else ref['reason'][:100]}", flush=True)

    def block(self, name):
        return [sys.executable, "-I", self.blocks / name]

    def priva(self, *words):
        return [self.binary, *words]

    # ------------------------------------------------------------- CDS FASTA
    def cds(self):
        good = b">ENST00000000233.10 cds chromosome:GRCh38:7:127588411:127591700:1 gene:ENSG1\nATGGCC\nTTGA\n>ENST2 cds\nATG\n"
        cases = [
            ("cds_good_gzip", "x.fa.gz", good, "hg38", True),
            ("cds_good_plain", "x.fa", good, "GRCh38", False),
            ("cds_conflict", "x.fa.gz", good, "hg19", True),
            ("cds_two_assemblies", "x.fa", b">a GRCh37 GRCh38\nA\n", "GRCh38", False),
            ("cds_grch37", "x.fa", b">a chromosome:GRCh37:1\nA\n", "GRCh37", False),
            ("cds_empty_sequence", "x.fa", b">a\n>b\nA\n", "hg38", False),
            ("cds_empty_header", "x.fa", b">   \nA\n", "hg38", False),
            ("cds_bare_header", "x.fa", b">\nA\n", "hg38", False),
            ("cds_sequence_first", "x.fa", b"ATG\n>a\nA\n", "hg38", False),
            ("cds_last_empty", "x.fa", b">a\nA\n>b\n", "hg38", False),
            ("cds_empty_file", "x.fa", b"", "hg38", False),
            ("cds_blank_only", "x.fa", b"\n \t\n\n", "hg38", False),
            ("cds_crlf", "x.fa", b">a GRCh38\r\nAT\r\n\r\n>b\r\nG\r\n", "hg38", False),
            ("cds_lone_cr", "x.fa", b">a\rATG\r>b\rC", "hg38", False),
            ("cds_word_boundary", "x.fa", b">a GRCh380 xGRCh37 GRCh37_b\nA\n", "hg38", False),
            ("cds_unicode_word", "x.fa", ">a αGRCh37 GRCh38é\nA\n".encode(), "hg38", False),
            ("cds_unicode_space", "x.fa", ">a　\nA \n \n".encode(), "hg38", False),
            ("cds_blank_between", "x.fa", b">a\n\n\nA\n\n>b\nC\n", "hg38", False),
            ("cds_hg19_name", "x.fa.gz", b">a GRCh37\nA\n", "hg19", True),
        ]
        for name, filename, text, assembly, compressed in cases:
            def build(root, filename=filename, text=text, compressed=compressed):
                path = root / filename
                if compressed:
                    gz(path, text)
                else:
                    path.write_bytes(text)
                return path
            self.case(name, build, lambda p, a=assembly: self.block("cds_fasta.py") + [p, a],
                      lambda p, a=assembly: self.priva("cds-fasta", "validate", "--fasta", p, "--assembly", a))
        # Python stops with a traceback here (KeyError, UnicodeDecodeError).
        self.case("cds_unknown_assembly", lambda r: (r / "x.fa").write_bytes(good) and r / "x.fa",
                  lambda p: self.block("cds_fasta.py") + [p, "GRCh39"],
                  lambda p: self.priva("cds-fasta", "validate", "--fasta", p, "--assembly", "GRCh39"), message=False)
        self.case("cds_not_utf8", lambda r: (r / "x.fa").write_bytes(b">a\n\xff\xfe\n") and r / "x.fa",
                  lambda p: self.block("cds_fasta.py") + [p, "hg38"],
                  lambda p: self.priva("cds-fasta", "validate", "--fasta", p, "--assembly", "hg38"), message=False)
        multi = lambda r: (r / "m.fa.gz").write_bytes(gzip.compress(b">a GRCh38\n") + gzip.compress(b"ACGT\n")) and r / "m.fa.gz"
        self.case("cds_multi_member_gzip", multi, lambda p: self.block("cds_fasta.py") + [p, "hg38"],
                  lambda p: self.priva("cds-fasta", "validate", "--fasta", p, "--assembly", "hg38"))

    # ---------------------------------------------------------------- LOFTEE
    def loftee(self):
        def database(path: Path, table=True):
            connection = sqlite3.connect(path)
            if table:
                connection.execute("CREATE TABLE phylocsf_data (chrom TEXT, pos INT)")
                connection.executemany("INSERT INTO phylocsf_data VALUES (?, ?)", [("1", i) for i in range(500)])
            connection.commit()
            connection.close()
            return path

        hg38 = "1\t248956422\t3\t60\t61\n2\t242193529\t253105752\t60\t61\n"
        def build_with(fai=hg38, first=b">1\nACGT\n", conservation="table"):
            def build(root):
                ancestor = gz(root / "human_ancestor.fa.gz", first)
                Path(str(ancestor) + ".fai").write_bytes(fai.encode() if isinstance(fai, str) else fai)
                sql = root / "loftee.sql"
                if conservation == "table":
                    database(sql)
                elif conservation == "empty":
                    database(sql, table=False)
                elif conservation == "text":
                    sql.write_bytes(b"This is not an SQLite database, only text of some length." * 4)
                elif conservation == "truncated":
                    database(sql)
                    data = sql.read_bytes()
                    sql.write_bytes(data[: len(data) // 2])
                elif conservation == "link":
                    database(root / "real.sql")
                    sql.symlink_to(root / "real.sql")
                return ancestor, sql
            return build
        cases = [
            ("loftee_good_hg38", build_with(), "hg38"),
            ("loftee_good_grch38", build_with(), "GRCh38"),
            ("loftee_good_hg19", build_with(fai="chr1\t249250621\t3\t60\t61\n"), "hg19"),
            ("loftee_wrong_assembly", build_with(), "GRCh37"),
            ("loftee_other_assembly_name", build_with(), "anything"),
            ("loftee_fai_four_fields", build_with(fai="1\t248956422\t3\t60\n"), "hg38"),
            ("loftee_fai_text", build_with(fai="1\t248956422\tx\t60\t61\n"), "hg38"),
            ("loftee_fai_zero_length", build_with(fai="1\t0\t3\t60\t61\n"), "hg38"),
            ("loftee_fai_conflict", build_with(fai=hg38 + "chr1\t5\t3\t60\t61\n"), "hg38"),
            ("loftee_fai_same_twice", build_with(fai=hg38 + "chr1\t248956422\t3\t60\t61\n"), "hg38"),
            ("loftee_fai_blank_line", build_with(fai=hg38 + "\n"), "hg38"),
            ("loftee_fai_crlf", build_with(fai=hg38.replace("\n", "\r\n")), "hg38"),
            ("loftee_fai_extra_fields", build_with(fai="1\t248956422\t3\t60\t61\tx\ty\n"), "hg38"),
            ("loftee_fai_no_chr1", build_with(fai="2\t242193529\t3\t60\t61\n"), "hg38"),
            ("loftee_not_fasta", build_with(first=b"ACGT\n"), "hg38"),
            ("loftee_empty_ancestor", build_with(first=b""), "hg38"),
            ("loftee_no_tables", build_with(conservation="empty"), "hg38"),
            ("loftee_text_database", build_with(conservation="text"), "hg38"),
            ("loftee_truncated_database", build_with(conservation="truncated"), "hg38"),
            ("loftee_missing_database", build_with(conservation="missing"), "hg38"),
            ("loftee_linked_database", build_with(conservation="link"), "hg38"),
        ]
        for name, build, assembly in cases:
            self.case(name, build,
                      lambda c, a=assembly: self.block("loftee_resources.py") + [c[0], c[1], a],
                      lambda c, a=assembly: self.priva("loftee-resources", "validate", "--ancestor", c[0],
                                                       "--conservation", c[1], "--assembly", a))

    # ---------------------------------------------------- mechanism evidence
    COLUMNS = ["gene_symbol", "source", "source_record_id", "source_condition_id", "mondo_id", "disease_scope",
               "priva_scope", "scope_review_status", "disease_label", "inheritance", "penetrance_raw",
               "penetrance_hpo_ids", "normalized_penetrance", "patho_mode_raw", "normalized_mechanisms",
               "mechanism_confidence", "disease_confidence", "pmids"]

    def evidence(self):
        def row(source="G2P_DDG2P", condition="OMIM:1", mondo="", penetrance="", raw="", mechanisms="", **extra):
            values = dict(gene_symbol="GENE", source=source, source_condition_id=condition, mondo_id=mondo,
                          normalized_penetrance=penetrance, patho_mode_raw=raw, normalized_mechanisms=mechanisms)
            values.update(extra)
            return values
        g2p_lof = row(raw="loss of function", mechanisms="LOF")
        g2p_pen = row(penetrance="complete")
        orpha = row(source="Orphadata", condition="", mondo="MONDO:1", mechanisms="GOF")
        good = [g2p_lof, g2p_pen, orpha]

        def table(rows, columns=None, newline="\n", blank=False, raw_text=None):
            def build(root):
                path = root / "evidence.tsv"
                if raw_text is not None:
                    path.write_bytes(raw_text)
                    return path
                names = columns or self.COLUMNS
                lines = ["\t".join(names)]
                for values in rows:
                    lines.append("\t".join(values.get(name, "") for name in names))
                    if blank:
                        lines.append("")
                path.write_bytes((newline.join(lines) + newline).encode())
                return path
            return build
        header = "\t".join(self.COLUMNS)
        cases = [
            ("evidence_good", table(good)),
            ("evidence_blank_lines", table(good, blank=True)),
            ("evidence_crlf", table(good, newline="\r\n")),
            ("evidence_reordered_extra", table(good, columns=list(reversed(self.COLUMNS)) + ["extra"])),
            ("evidence_missing_columns", table(good, columns=[c for c in self.COLUMNS if c not in ("pmids", "inheritance")])),
            ("evidence_header_only", table([])),
            ("evidence_unexpected_sources", table(good + [row(source="ClinGen"), row(source="")])),
            ("evidence_no_g2p", table([orpha])),
            ("evidence_no_lof", table([g2p_pen, orpha])),
            ("evidence_lof_by_raw_case", table([row(raw="Loss Of Function"), g2p_pen, orpha])),
            ("evidence_lof_by_mechanism_case", table([row(mechanisms="lof"), g2p_pen, orpha])),
            ("evidence_no_penetrance", table([g2p_lof, orpha])),
            ("evidence_penetrance_spaces", table([g2p_lof, row(penetrance="  InComplete \t"), orpha])),
            ("evidence_bad_penetrance", table([g2p_lof, row(penetrance="Partial"), orpha])),
            ("evidence_bad_penetrance_quote", table([g2p_lof, row(penetrance="it's"), orpha])),
            ("evidence_penetrance_no_condition", table([g2p_lof, row(condition="", penetrance="complete"), orpha])),
            ("evidence_penetrance_mondo_only", table([g2p_lof, row(condition="", mondo="MONDO:2", penetrance="complete"), orpha])),
            ("evidence_no_orphadata", table([g2p_lof, g2p_pen])),
            ("evidence_no_orphadata_mondo", table([g2p_lof, g2p_pen, row(source="Orphadata", condition="")])),
            ("evidence_quoted_tab", table(good, raw_text=(header + "\n" + "\t".join(
                ['"GE\tNE"' if c == "gene_symbol" else g2p_lof.get(c, "") for c in self.COLUMNS]) + "\n" + "\n".join(
                "\t".join(v.get(c, "") for c in self.COLUMNS) for v in (g2p_pen, orpha)) + "\n").encode())),
            ("evidence_short_orphadata_row", table(good, raw_text=(header + "\n" + "\n".join(
                "\t".join(v.get(c, "") for c in self.COLUMNS) for v in good) + "\nGENE\tOrphadata\n").encode())),
            ("evidence_empty_file", table([], raw_text=b"")),
        ]
        for name, build in cases:
            self.case(name, build, lambda p: self.block("mechanism_evidence.py") + [p],
                      lambda p: self.priva("mechanism-evidence", "validate", "--evidence", p))
        # A G2P row shorter than its header: Python stops with AttributeError.
        short = table(good, raw_text=(header + "\n" + "\n".join(
            "\t".join(v.get(c, "") for c in self.COLUMNS) for v in good) + "\nGENE\tG2P_DDG2P\n").encode())
        self.case("evidence_short_g2p_row", short, lambda p: self.block("mechanism_evidence.py") + [p],
                  lambda p: self.priva("mechanism-evidence", "validate", "--evidence", p), message=False)

    # ----------------------------------------------------------- MONDO scope
    def mondo(self):
        registry_columns = ["disease_id", "mondo_id", "disease_scope", "priva_scope", "scope_evidence", "scope_review_status", "note"]
        assertion_columns = ["gene_symbol", "disease_id", "hpo_id", "frequency", "evidence", "reference", "mondo_id",
                             "disease_scope", "priva_scope", "scope_evidence", "scope_review_status"]

        def tables(registry_rows=2, assertion_rows=3, registry_gz=False, assertions_gz=True,
                   registry_cols=registry_columns, assertion_cols=assertion_columns, newline="\n", blank=False):
            def build(root):
                def write(name, columns, count, compressed):
                    lines = ["\t".join(columns)]
                    for index in range(count):
                        lines.append("\t".join(f"{c}{index}" for c in columns))
                        if blank:
                            lines.append("")
                    data = (newline.join(lines) + newline).encode()
                    path = root / (name + (".tsv.gz" if compressed else ".tsv"))
                    path.write_bytes(gzip.compress(data) if compressed else data)
                    return path
                return (write("registry", registry_cols, registry_rows, registry_gz),
                        write("assertions", assertion_cols, assertion_rows, assertions_gz))
            return build
        cases = [
            ("mondo_good", tables()),
            ("mondo_both_gzip", tables(registry_gz=True)),
            ("mondo_both_plain", tables(assertions_gz=False)),
            ("mondo_crlf_blank", tables(newline="\r\n", blank=True)),
            ("mondo_registry_missing", tables(registry_cols=["disease_id", "note"])),
            ("mondo_assertions_missing", tables(assertion_cols=assertion_columns[:-3])),
            ("mondo_registry_empty", tables(registry_rows=0)),
            ("mondo_assertions_empty", tables(assertion_rows=0)),
            ("mondo_no_header", tables(registry_cols=[], registry_rows=0)),
        ]
        for name, build in cases:
            self.case(name, build, lambda c: self.block("mondo_scope.py") + [c[0], c[1]],
                      lambda c: self.priva("mondo-scope", "validate", "--registry", c[0], "--hpo-assertions", c[1]))

    # ------------------------------------------------------- SpliceAI2 check
    @staticmethod
    def manifest(**changes):
        value = dict(status="downloads_complete", confirmed_by="owner", assembly="GRCh38", version="2.0",
                     source_revision=REVISION, files=[{"path": "/data/spliceai2/chr1.tsv.gz", "sha256": "0" * 64}])
        value.update(changes)
        return value

    def spliceai2_check(self):
        def setup(config="assembly: hg38\nthreads: 4\n", manifest=None, output="absent", manifest_name="manifest.json",
                  link_manifest=False, raw_manifest=None):
            def build(root):
                (root / "config.yaml").write_text(config)
                folder = root / "source"
                folder.mkdir(parents=True, exist_ok=True)
                target = root / manifest_name
                target.parent.mkdir(parents=True, exist_ok=True)
                target.write_text(raw_manifest if raw_manifest is not None else json.dumps(manifest or self.manifest()))
                if link_manifest:
                    (root / "linked.json").symlink_to(target)
                    target = root / "linked.json"
                out = root / "out"
                if output in ("empty", "building", "other", "ready"):
                    out.mkdir()
                if output == "building":
                    (out / ".building-123").mkdir()
                elif output == "other":
                    (out / "chr1.vcf.gz").write_text("x")
                elif output == "ready":
                    (out / "chr1.vcf.gz").write_text("x")
                    (out / "adapter.ready.json").write_text("{}")
                elif output == "file":
                    out.write_text("not a folder")
                return root / "config.yaml", target, out
            return build
        cases = [
            ("check_good", setup(), None),
            ("check_grch38", setup(config="assembly: GRCh38\n"), None),
            ("check_hg19", setup(config="assembly: hg19\n"), None),
            ("check_no_assembly", setup(config="threads: 4\n"), None),
            ("check_config_list", setup(config="- a\n- b\n"), None),
            ("check_config_empty", setup(config=""), None),
            ("check_relative_manifest", setup(), "relative"),
            ("check_partial_manifest", setup(manifest_name="dl.PARTIAL/manifest.json"), None),
            ("check_linked_manifest", setup(link_manifest=True), None),
            ("check_parent_component", setup(), "dotdot"),
            ("check_dot_component", setup(), "dot"),
            ("check_trailing_slash_output", setup(output="empty"), "slash"),
            ("check_status", setup(manifest=self.manifest(status="downloading")), None),
            ("check_blank_owner", setup(manifest=self.manifest(confirmed_by="  ")), None),
            ("check_owner_number", setup(manifest=self.manifest(confirmed_by=7)), None),
            ("check_revision", setup(manifest=self.manifest(source_revision="0" * 40)), None),
            ("check_version_number", setup(manifest=self.manifest(version=2.0)), None),
            ("check_partial_source", setup(manifest=self.manifest(files=[{"nested": [{"path": "/x/a.partial"}]}])), None),
            ("check_path_number", setup(manifest=self.manifest(files=[{"path": 5}])), None),
            ("check_not_object", setup(raw_manifest="[1, 2]"), None),
            ("check_output_empty", setup(output="empty"), None),
            ("check_output_building", setup(output="building"), None),
            ("check_output_other", setup(output="other"), None),
            ("check_output_ready", setup(output="ready"), None),
            ("check_output_file", setup(output="file"), None),
        ]
        for name, build, variant in cases:
            def argv(c, variant=variant):
                config, manifest, output = c
                if variant == "relative":
                    manifest = os.path.relpath(manifest)
                elif variant == "dotdot":
                    manifest = f"{manifest.parent}/../{manifest.parent.name}/{manifest.name}"
                elif variant == "dot":
                    manifest = f"{manifest.parent}/./{manifest.name}"
                elif variant == "slash":
                    output = f"{output}/"
                return config, manifest, output
            self.case(name, build, lambda c, argv=argv: self.block("spliceai2_check.py") + list(argv(c)),
                      lambda c, argv=argv: self.priva("spliceai2-install", "check", *sum(
                          ([f"--{k}", v] for k, v in zip(("config", "manifest", "output-dir"), argv(c))), [])))
        self.case("check_manifest_not_json", setup(raw_manifest="{broken"),
                  lambda c: self.block("spliceai2_check.py") + list(c),
                  lambda c: self.priva("spliceai2-install", "check", "--config", c[0], "--manifest", c[1], "--output-dir", c[2]),
                  message=False)

    # ---------------------------------------------------- SpliceAI2 register
    def spliceai2_register(self):
        config_text = ("# PriVA configuration\nassembly: hg38\nspliceai_snv_prescore: /old/snv.vcf.gz  # v1\n"
                       "threads: 4\nnested:\n  spliceai_indel_prescore: kept\nspliceai_indel_prescore: \"/old/indel.vcf.gz\"\n")

        def setup(config=config_text, receipt_change=None, output_change=None, folder="out", tamper=None):
            def build(root):
                out = root / folder
                out.mkdir(parents=True)
                (root / "config.yaml").write_text(config) if isinstance(config, str) else (root / "config.yaml").write_bytes(config)
                manifest = root / "manifest.json"
                manifest.write_text(json.dumps(self.manifest()))
                outputs = {}
                for role in ("snv", "indel"):
                    vcf = out / f"spliceai2.{role}.vcf.gz"
                    vcf.write_bytes(gzip.compress(f"##fileformat=VCFv4.2 {role}\n".encode()))
                    index = Path(str(vcf) + ".tbi")
                    index.write_bytes(b"TBI\x01" + role.encode())
                    outputs[role] = dict(path=str(vcf), sha256=sha(vcf.read_bytes()), bytes=vcf.stat().st_size,
                                         index_path=str(index), index_sha256=sha(index.read_bytes()))
                receipt = dict(schema_version=1, resource="SpliceAI2", version="2.0", assembly="GRCh38",
                               status="prepared_annotation_only", scientific_qualified=False, coverage="complete_published",
                               source_revision=REVISION, source_manifest_sha256=sha(manifest.read_bytes()), outputs=outputs)
                if receipt_change:
                    receipt_change(receipt, out)
                if output_change:
                    output_change(out, outputs)
                path = out / "adapter.ready.json"
                path.write_text(json.dumps(receipt, indent=2))
                context = dict(config=root / "config.yaml", manifest=manifest, out=out, receipt=path,
                               config_sha=sha((root / "config.yaml").read_bytes()), manifest_sha=sha(manifest.read_bytes()))
                if tamper:
                    tamper(context)
                return context
            return build

        def args(c):
            return [c["config"], c["manifest"], c["out"], c["receipt"], c["config_sha"], c["manifest_sha"]]

        def native(c):
            names = ("config", "manifest", "output-dir", "receipt", "config-sha", "manifest-sha")
            return self.priva("spliceai2-install", "register", *sum(([f"--{n}", v] for n, v in zip(names, args(c))), []))

        def set_key(key, value):
            return lambda receipt, out: receipt.update({key: value})

        def outputs_change(change):
            return lambda receipt, out: change(receipt["outputs"], out)

        def link_folder(root_name):
            def tamper(c):
                linked = c["out"].parent / "linked"
                linked.symlink_to(c["out"])
                c["out"] = linked
                c["receipt"] = linked / "adapter.ready.json"
            return tamper
        cases = [
            ("register_replace", setup(), {}),
            ("register_append", setup(config="assembly: hg38\nthreads: 4"), {}),
            ("register_crlf", setup(config=config_text.replace("\n", "\r\n")), {}),
            ("register_unicode_folder", setup(folder="café/out"), {}),
            ("register_space_before_colon", setup(config="assembly: hg38\nspliceai_snv_prescore : x\n"), {}),
            ("register_duplicate_key", setup(config="assembly: hg38\nspliceai_snv_prescore: a\nspliceai_snv_prescore: b\n"), {}),
            ("register_config_changed", setup(tamper=lambda c: c["config"].write_text(config_text + "x: 1\n")), {}),
            ("register_manifest_changed", setup(tamper=lambda c: c["manifest"].write_text("{}")), {}),
            ("register_schema_bool", setup(receipt_change=set_key("schema_version", True)), {}),
            ("register_schema_float", setup(receipt_change=set_key("schema_version", 1.0)), {}),
            ("register_qualified_zero", setup(receipt_change=set_key("scientific_qualified", 0)), {}),
            ("register_status", setup(receipt_change=set_key("status", "prepared")), {}),
            ("register_manifest_sha", setup(receipt_change=set_key("source_manifest_sha256", "0" * 64)), {}),
            ("register_not_object", setup(tamper=lambda c: c["receipt"].write_text("[]")), {}),
            ("register_outputs_one", setup(receipt_change=outputs_change(lambda o, out: o.pop("indel"))), {}),
            ("register_outputs_extra", setup(receipt_change=outputs_change(lambda o, out: o.update(other={}))), {}),
            ("register_output_list", setup(receipt_change=outputs_change(lambda o, out: o.update(snv=[]))), {}),
            ("register_checksum_upper", setup(receipt_change=outputs_change(lambda o, out: o["snv"].update(sha256=o["snv"]["sha256"].upper()))), {}),
            ("register_wrong_sha", setup(receipt_change=outputs_change(lambda o, out: o["indel"].update(sha256="0" * 64))), {}),
            ("register_wrong_index_sha", setup(receipt_change=outputs_change(lambda o, out: o["snv"].update(index_sha256="1" * 64))), {}),
            ("register_wrong_bytes", setup(receipt_change=outputs_change(lambda o, out: o["snv"].update(bytes=1))), {}),
            ("register_bytes_text", setup(receipt_change=outputs_change(lambda o, out: o["snv"].update(bytes=str(o["snv"]["bytes"])))), {}),
            ("register_index_name", setup(receipt_change=outputs_change(lambda o, out: o["snv"].update(index_path=o["snv"]["path"] + ".csi"))), {}),
            ("register_same_file_twice", setup(receipt_change=outputs_change(lambda o, out: o.update(indel=dict(o["snv"])))), {}),
            ("register_outside", setup(receipt_change=outputs_change(lambda o, out: o["snv"].update(path=str(out.parent / "x.vcf.gz")))), {}),
            ("register_relative_path", setup(receipt_change=outputs_change(lambda o, out: o["snv"].update(path="out/x.vcf.gz"))), {}),
            ("register_dot_path", setup(receipt_change=outputs_change(lambda o, out: o["snv"].update(
                path=str(out / "." / Path(o["snv"]["path"]).name), index_path=str(out / "." / Path(o["snv"]["path"]).name) + ".tbi"))), {}),
            ("register_building", setup(tamper=lambda c: (c["out"] / ".building-1").mkdir()), {}),
            ("register_linked_folder", setup(tamper=link_folder("out")), {}),
            ("register_partial_folder", setup(folder="x.partial/out"), {}),
        ]
        for name, build, _ in cases:
            self.case(name, build, lambda c: self.block("spliceai2_register.py") + args(c), native,
                      written=lambda c: c["config"])
        vcf_name = setup(output_change=lambda out, outputs: None,
                         receipt_change=outputs_change(lambda o, out: None))
        def plain_vcf(receipt, out):
            old = Path(receipt["outputs"]["snv"]["path"])
            new = out / "spliceai2.snv.vcf"
            old.rename(new)
            Path(str(old) + ".tbi").rename(str(new) + ".tbi")
            receipt["outputs"]["snv"].update(path=str(new), index_path=str(new) + ".tbi")
        self.case("register_vcf_name", setup(receipt_change=plain_vcf),
                  lambda c: self.block("spliceai2_register.py") + args(c), native, written=lambda c: c["config"])
        self.case("register_lock_held", setup(), lambda c: self.block("spliceai2_register.py") + args(c), native,
                  written=lambda c: c["config"], hold_lock=lambda c: c["out"] / ".install.lock")

    # -------------------------------------------------------------- VEP pins
    def pins(self):
        def records(*items):
            return list(items)
        base = [
            {"name": "zlib", "version": "1.3.1", "build_string": "hb9d3cd8_2", "channel": "conda-forge", "platform": "linux-64"},
            {"name": "ensembl-vep", "version": "113.2", "build_string": "pl5321h2a3209d_0", "channel": "bioconda"},
            {"name": "diptest", "version": "0.8.1", "build_string": "pypi_0", "channel": "pypi"},
            {"name": "bcftools", "version": "1.21", "build_string": "h3a4d415_1", "channel": "bioconda"},
            {"name": "_libgcc_mutex", "version": "0.1", "build_string": "conda_forge", "channel": "conda-forge"},
        ]
        cases = [
            ("pins_micromamba_list", base),
            ("pins_conda_object", {"packages": base}),
            ("pins_build_fallback", [{"name": "a", "version": "1", "build_string": "", "build": "h0"}]),
            ("pins_platform_pypi", base + [{"name": "p", "version": "1", "build": "0", "platform": "pypi"}]),
            ("pins_bad_version", base + [{"name": "a", "version": "1 2", "build": "0"}]),
            ("pins_number_version", base + [{"name": "a", "version": 1, "build": "0"}]),
            ("pins_no_build", base + [{"name": "a", "version": "1"}]),
            ("pins_duplicate", base + [dict(base[0])]),
            ("pins_only_vep", [base[1], base[2]]),
            ("pins_empty", []),
            ("pins_object_without_packages", {"other": base}),
            ("pins_not_records", [1, 2]),
            ("pins_unicode_name", base + [{"name": "café", "version": "1", "build": "0"}]),
        ]
        for name, value in cases:
            def build(root, value=value):
                path = root / "pins.txt"
                path.write_text(json.dumps(value))
                return path
            self.case(name, build, lambda p: self.block("vep_release_pins.py") + [p],
                      lambda p: self.priva("vep-release", "pins", "--packages", p), written=lambda p: p)
        for name, text in [("pins_missing_name", json.dumps([{"version": "1", "build": "0"}])),
                           ("pins_missing_version", json.dumps([{"name": "a", "build": "0"}]))]:
            self.case(name, lambda r, text=text: (r / "pins.txt").write_text(text) and r / "pins.txt",
                      lambda p: self.block("vep_release_pins.py") + [p],
                      lambda p: self.priva("vep-release", "pins", "--packages", p), written=lambda p: p)
        self.case("pins_not_json", lambda r: (r / "pins.txt").write_text("{") and r / "pins.txt",
                  lambda p: self.block("vep_release_pins.py") + [p],
                  lambda p: self.priva("vep-release", "pins", "--packages", p), written=lambda p: p, message=False)


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = parser.add_subparsers(dest="command", required=True)
    check = sub.add_parser("check")
    check.add_argument("--repo", type=Path, required=True)
    check.add_argument("--binary", type=Path, required=True)
    check.add_argument("--work-root", type=Path, required=True)
    args = parser.parse_args()
    args.work_root.mkdir(parents=True, exist_ok=False)
    compare = Compare(args.repo.resolve(), args.binary.resolve(), args.work_root.resolve())
    os.chdir(args.work_root)
    for part in (compare.cds, compare.loftee, compare.evidence, compare.mondo,
                 compare.spliceai2_check, compare.spliceai2_register, compare.pins):
        part()
    mismatches = [row for row in compare.rows if not row["passed"]]
    report = dict(cases=len(compare.rows), passed=len(compare.rows) - len(mismatches),
                  messages_checked=sum(row["message_checked"] for row in compare.rows),
                  successes=sum(row["reference"]["success"] for row in compare.rows), mismatches=mismatches)
    (args.work_root / "report.json").write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")
    print(json.dumps({key: value for key, value in report.items() if key != "mismatches"} | {"mismatch_count": len(mismatches)}))
    raise SystemExit(1 if mismatches else 0)


if __name__ == "__main__":
    main()
