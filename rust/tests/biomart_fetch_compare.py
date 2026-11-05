#!/usr/bin/env python3
"""Compare `priva biomart fetch` with biomart_nm_to_enst.py `fetch`.

  python biomart_fetch_compare.py check --binary PRIVA --repo CHECKOUT \
      --raw-dumps DIR --work-root NEW_DIR [--jobs N] [--only CASE ...] [--skip-full]

Neither program is changed for the comparison and neither reaches the network.
Each arm of each case gets its own local publisher: an HTTPS proxy that answers
CONNECT for the real Ensembl hosts, ends the TLS session with a local test
authority (trusted through SSL_CERT_FILE and CURL_CA_BUNDLE) and serves real
bytes, so the pinned URLs, sizes and SHA-256 values stay in force. The service
responses are rebuilt from the bundled tables and the registry and dataset
text in their receipts; the FTP files are the retained dumps in --raw-dumps.
A case's plan makes single responses fail: a status, a cut body, a reset, a
redirect, an ignored Range.

Exit status, failure message, stdout, the requests the publisher received and
every file left in the arm are compared, file bytes by SHA-256. In JSON
receipts only the run-dependent text is masked: arm path, timestamps, tempfile
tokens and the producer fields, which are checked separately for each program.
After a success, each program verifies the other's output, `identity` must
agree, and a second fetch must reuse the table without a request.
"""
from __future__ import annotations

import argparse
import concurrent.futures
import fcntl
import hashlib
import http.server
import json
import os
import re
import shutil
import socket
import socketserver
import ssl
import struct
import subprocess
import sys
import threading
import time
from pathlib import Path

G37 = "ensembl-grch37-mart116"
G38 = "ensembl-grch38-mart113"
FTP = "ensembl-grch38-mart116-ftp"
BASENAME = {G37: "ensembl_GRCh37_mart116_nm_to_enst.tsv",
            G38: "ensembl_GRCh38_mart113_nm_to_enst.tsv",
            FTP: "ensembl_GRCh38_mart116_ftp_candidate.tsv"}
ASSEMBLY = {G37: "GRCh37", G38: "GRCh38", FTP: "GRCh38"}
HOSTS = ("grch37.ensembl.org", "oct2024.archive.ensembl.org", "ftp.ensembl.org",
         "elsewhere.example.org")
NATIVE_PRODUCER = "priva biomart fetch"
FTP_FILES = ("CHECKSUMS", "ensembl_mart_116.sql.gz", "meta_conf__dataset__main.txt.gz",
             "meta_conf__xml__dm.txt.gz", "hsapiens_gene_ensembl__transcript__main.txt.gz",
             "hsapiens_gene_ensembl__translation__main.txt.gz", "hsapiens_gene_ensembl__ox_hgnc__dm.txt.gz",
             "hsapiens_gene_ensembl__ox_uniprotswissprot__dm.txt.gz",
             "hsapiens_gene_ensembl__ox_refseq_mrna__dm.txt.gz")
FULL = ("f_full_cold", "f_full_resume", "f_full_canonical", "f_full_retained_dumps")


def sha256_file(path):
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def make_authority(folder):
    """A test authority and one server certificate for the publisher hosts."""
    folder.mkdir()
    run = lambda *command: subprocess.run(command, cwd=folder, check=True, capture_output=True)
    run("openssl", "req", "-x509", "-newkey", "rsa:2048", "-nodes", "-keyout", "ca.key", "-out", "ca.pem",
        "-days", "3", "-subj", "/CN=local test authority",
        "-addext", "basicConstraints=critical,CA:TRUE", "-addext", "keyUsage=critical,keyCertSign")
    run("openssl", "req", "-newkey", "rsa:2048", "-nodes", "-keyout", "leaf.key", "-out", "leaf.csr",
        "-subj", "/CN=local publisher")
    (folder / "leaf.ext").write_text(
        "subjectAltName=" + ",".join("DNS:" + host for host in HOSTS) + "\n"
        "extendedKeyUsage=serverAuth\nkeyUsage=critical,digitalSignature,keyEncipherment\n")
    run("openssl", "x509", "-req", "-in", "leaf.csr", "-CA", "ca.pem", "-CAkey", "ca.key", "-CAcreateserial",
        "-out", "leaf.pem", "-days", "3", "-extfile", "leaf.ext")
    context = ssl.SSLContext(ssl.PROTOCOL_TLS_SERVER)
    context.load_cert_chain(folder / "leaf.pem", folder / "leaf.key")
    return context, folder / "ca.pem"


class SocketWriter:
    """The handler's output: every byte goes out, however TLS records split it."""
    closed = False

    def __init__(self, connection):
        self.connection = connection

    def write(self, data):
        self.connection.sendall(data)
        return len(data)

    def flush(self):
        pass

    def close(self):
        pass


class Publisher(socketserver.ThreadingTCPServer):
    allow_reuse_address = True
    daemon_threads = True

    def handle_error(self, request, client_address):
        pass  # A planned reset ends its handler in the middle of a reply.

    def __init__(self, tls, content, plan):
        super().__init__(("127.0.0.1", 0), Reply)
        self.tls = tls
        self.content = content
        self.plan = {key: list(value) for key, value in (plan or {}).items()}
        self.requests = []
        self.lock = threading.Lock()

    def next_action(self, kind):
        with self.lock:
            queue = self.plan.get(kind)
            return queue.pop(0) if queue else "ok"

    def __enter__(self):
        threading.Thread(target=self.serve_forever, daemon=True).start()
        return self

    def __exit__(self, *_):
        self.shutdown()
        self.server_close()


class Reply(http.server.BaseHTTPRequestHandler):
    protocol_version = "HTTP/1.1"
    timeout = 120

    def log_message(self, *_):
        pass

    def finish(self):
        try:
            super().finish()
            self.connection.close()
        except (OSError, ValueError):
            pass

    def do_CONNECT(self):
        self.send_response(200, "Connection established")
        self.end_headers()
        self.host = self.path.split(":")[0]
        try:
            self.connection = self.server.tls.wrap_socket(self.connection, server_side=True)
        except (OSError, ssl.SSLError):
            self.close_connection = True
            return
        self.rfile = self.connection.makefile("rb", -1)
        self.wfile = SocketWriter(self.connection)
        self.close_connection = False

    def kind(self):
        name = self.path.rsplit("/", 1)[-1]
        if name in FTP_FILES:
            return name
        for marker, kind in (("type=registry", "registry"), ("type=datasets", "datasets"), ("query=", "export")):
            if marker in self.path:
                return kind
        return "unknown"

    def do_GET(self):
        kind = self.kind()
        self.server.requests.append(dict(
            host=self.host, path=self.path, range=self.headers.get("Range"),
            accept_encoding=self.headers.get("Accept-Encoding"), user_agent=self.headers.get("User-Agent")))
        action = self.server.next_action(kind)
        body = self.server.content(kind)
        if body is None:
            action = {"status": 404}
        elif isinstance(action, dict) and "body" in action:
            body = action["body"]
        self.close_connection = True
        self.respond(action, body)

    def head(self, status, reason, *fields):
        self.send_response_only(status, reason)
        for name, value in fields:
            self.send_header(name, value)
        self.end_headers()

    def close_clean(self):
        """End the TLS session properly: the peer reads end of data."""
        try:
            self.connection.settimeout(2)
            self.connection = self.wfile.connection = self.connection.unwrap()
        except (OSError, ValueError):
            pass

    def close_reset(self):
        self.connection.setsockopt(socket.SOL_SOCKET, socket.SO_LINGER, struct.pack("ii", 1, 0))
        self.connection.close()

    def respond(self, action, body):
        action = {action: True} if isinstance(action, str) else action
        wanted = re.fullmatch(r"bytes=(\d+)-", self.headers.get("Range") or "")
        if "status" in action:
            reason = {404: "Not Found", 500: "Internal Server Error", 503: "Service Unavailable"}[action["status"]]
            self.head(action["status"], reason, ("Content-Length", "0"))
        elif "redirect" in action:
            self.head(302, "Found", ("Location", action["redirect"]), ("Content-Length", "0"))
        elif "ignore_range" in action or not wanted:
            self.send_body(action, body, 200, "OK", ())
        else:
            start = int(wanted[1])
            content_range = "bytes 0-1/2" if "bad_range" in action else f"bytes {start}-{len(body) - 1}/{len(body)}"
            self.send_body(action, body[start:], 206, "Partial Content", (("Content-Range", content_range),))

    def send_body(self, action, body, status, reason, fields):
        if "chunk_cut" in action:
            self.head(status, reason, ("Transfer-Encoding", "chunked"), *fields)
            piece = body[:action["chunk_cut"]]
            self.wfile.write(b"%x\r\n" % len(piece) + piece + b"\r\n")
            self.close_clean()
            return
        length = len(body) + action.get("extra_length", 0)
        if "no_length" in action:
            self.head(status, reason, *fields)
        else:
            self.head(status, reason, ("Content-Length", str(length)), *fields)
        cut = next((action[key] for key in ("cut", "abrupt", "reset") if key in action), None)
        self.wfile.write(body if cut is None else body[:cut])
        if "cut" in action:
            self.close_clean()
        elif "abrupt" in action:
            self.connection.shutdown(socket.SHUT_RDWR)
        elif "reset" in action:
            self.close_reset()


class Sources:
    """Real publisher bytes, read once."""

    def __init__(self, repo, raw_dumps):
        self.service = {}
        self.raw_dumps = raw_dumps
        for source in (G37, G38):
            table = repo / "data/hgnc" / BASENAME[source]
            receipt = json.loads(Path(str(table) + ".provenance.json").read_text())
            export = table.read_bytes() + b"[success]\n"
            assert hashlib.sha256(export).hexdigest() == receipt["raw_sha256"]
            self.service[source] = dict(
                registry=receipt["metadata_before"]["registry_xml"].encode(),
                datasets=receipt["metadata_before"]["datasets_tsv"].encode(),
                export=export)
        self.files = {}
        self.lock = threading.Lock()

    def content(self, source):
        def read(kind):
            if kind in FTP_FILES:
                with self.lock:
                    if kind not in self.files:
                        self.files[kind] = (self.raw_dumps / kind).read_bytes()
                return self.files[kind]
            return self.service.get(source, {}).get(kind)
        return read


def case(name, source=G37, **fields):
    return name, dict(source=source, **fields)


def copy_mapping(root, repo, source, name=None, receipt=True):
    folder = root / "out"
    folder.mkdir(exist_ok=True)
    table = repo / "data/hgnc" / BASENAME[source]
    target = folder / (name or BASENAME[source])
    shutil.copyfile(table, target)
    if receipt:
        shutil.copyfile(str(table) + ".provenance.json", str(target) + ".provenance.json")
    return target


def seed_dumps(root, raw_dumps, names):
    raw = root / "res/raw"
    raw.mkdir(parents=True, exist_ok=True)
    for name in names:
        (raw / name).symlink_to(raw_dumps / name)


def write_part(root, raw_dumps, name, size=1000):
    raw = root / "res/raw"
    raw.mkdir(parents=True, exist_ok=True)
    with (raw_dumps / name).open("rb") as handle:
        (raw / (name + ".part")).write_bytes(handle.read(size))


def changed_registry(xml):
    return re.sub(rb'(<MartURLLocation\b[^>]*\bname="ENSEMBL_MART_ENSEMBL")', rb'\1 extra="1"', xml)


def cases():
    """Each case: the source, the state left in the arm, and the publisher's plan.

    setup(root, context) prepares the arm. `transport` marks a failure whose
    text comes from the transport layer, so only its presence is compared.
    `volatile` marks a cut whose byte position depends on read buffering, so
    request offsets and retained partial bytes are not compared.
    """
    other = "https://elsewhere.example.org/biomart/martservice?type=registry&virtualSchema=default"
    same = "/biomart/martservice?type=registry&virtualSchema=default&moved=1"

    def existing(source, **options):
        return lambda root, context: copy_mapping(root, context.repo, source, **options)

    def link_part(root, _):
        (root / "out").mkdir()
        (root / "out" / (BASENAME[G37] + ".part")).symlink_to(root / "elsewhere")

    def raw_link_part(root, _):
        (root / "raw").mkdir()
        (root / "raw" / (BASENAME[G37] + ".raw.tsv.part")).symlink_to(root / "elsewhere")

    def receipt_only(source):
        def setup(root, context):
            target = copy_mapping(root, context.repo, source)
            target.unlink()
        return setup

    def raw_present(root, _):
        (root / "raw").mkdir()
        (root / "raw" / (BASENAME[G37] + ".raw.tsv")).write_text("published raw\n")

    def part(name, **options):
        return lambda root, context: write_part(root, context.raw_dumps, name, **options)

    def junk_raw(size):
        def setup(root, _):
            (root / "res/raw").mkdir(parents=True)
            (root / "res/raw/CHECKSUMS").write_bytes(b"x" * size)
        return setup

    def link_raw_part(root, _):
        (root / "res/raw").mkdir(parents=True)
        (root / "res/raw/CHECKSUMS.part").symlink_to(root / "elsewhere")

    def long_part(root, _):
        (root / "res/raw").mkdir(parents=True)
        (root / "res/raw/CHECKSUMS.part").write_bytes(b"x" * 3377161)

    def retained_dumps(root, context):
        seed_dumps(root, context.raw_dumps, FTP_FILES)

    def resume_seed(root, context):
        write_part(root, context.raw_dumps, FTP_FILES[4], 1_000_000)

    def bad_table(root, _):
        (root / "out").mkdir()
        (root / "out" / BASENAME[FTP]).write_text("not the table\n")

    def file_as_folder(root, _):
        (root / "res").write_text("file\n")

    not_found = {"status": 404}
    service_failed = {"status": 503}
    items = [
        # Clean real acquisitions, reuse and the bundled receipt.
        case("s_clean_grch37", G37),
        case("s_clean_grch38", G38),
        case("s_bundled_grch37", G37, output="source/data/hgnc/" + BASENAME[G37]),
        case("s_relative_paths", G37, relative=True),
        case("s_custom_name", G38, output="out/selected_mapping.tsv"),
        case("s_reuse", G37, setup=existing(G37)),
        case("s_reuse_wrong_source", G38, setup=existing(G37, name="selected.tsv"), output="out/selected.tsv"),
        case("s_reuse_table_without_receipt", G37, setup=existing(G37, receipt=False)),
        # Refusals before any request.
        case("s_name_conflict", G37, output="out/" + BASENAME[G38]),
        case("s_receipt_without_table", G37, setup=receipt_only(G37)),
        case("s_raw_exists", G37, setup=raw_present),
        case("s_raw_is_output", G37, raw_output="out/" + BASENAME[G37]),
        case("s_raw_in_checkout", G37, raw_output="source/raw/export.tsv"),
        case("s_resource_dir_refused", G37, extra=["--resource-dir", "res"]),
        case("s_lock_held", G37, hold_lock="out"),
        case("s_table_part_symlink", G37, setup=link_part),
        case("s_raw_part_symlink", G37, setup=raw_link_part),
        case("s_attempts_zero", G37, attempts=0),
        case("s_attempts_six", G37, attempts=6),
        case("s_timeout_zero", G37, timeout=0),
        case("s_no_raw_output", G37, raw_output=None, usage=True),
        case("s_unknown_source", "unknown-source", usage=True, output="out/x.tsv"),
        # Metadata requests.
        case("s_metadata_503_then_ok", G37, plan={"registry": [service_failed]}),
        case("s_metadata_503_exhausted", G37, plan={"datasets": [service_failed, service_failed]}),
        case("s_metadata_reset_then_ok", G37, plan={"registry": [{"reset": 10}]}, transport=True),
        case("s_metadata_redirect_away", G37, plan={"registry": [{"redirect": other}, "ok", {"redirect": other}, "ok"]}),
        case("s_metadata_redirect_same_host", G37, plan={"registry": [{"redirect": same}]}),
        case("s_metadata_bad_xml", G37, plan={"registry": [{"body": b"<broken"}]}),
        case("s_metadata_changed_source", G37, plan={"registry": [{"body": b"<x/>"}]}),
        case("s_metadata_not_utf8", G37, plan={"registry": [{"body": b"\xff<x/>"}, {"body": b"\xff<x/>"}]}),
        case("s_metadata_cut_body", G37, plan={"datasets": [{"cut": 20}]}),
        case("s_metadata_chunk_cut", G37, plan={"datasets": [{"chunk_cut": 20}]}, transport=True),
        # The export request.
        case("s_export_503_then_ok", G37, plan={"export": [service_failed]}),
        case("s_export_503_exhausted", G37, plan={"export": [service_failed, service_failed]}),
        case("s_export_404_one_attempt", G37, attempts=1, plan={"export": [not_found]}),
        case("s_export_reset_then_ok", G37, plan={"export": [{"reset": 3_000_000}]}, transport=True),
        case("s_export_cut_body", G37, plan={"export": [{"cut": 3_000_000}]}),
        case("s_export_abrupt_close", G37, plan={"export": [{"abrupt": 3_000_000}]}),
        case("s_export_chunk_cut_then_ok", G37, plan={"export": [{"chunk_cut": 2_000_000}]}, transport=True,
             volatile=True),
        case("s_export_redirect_away", G37, plan={"export": [{"redirect": other}]}),
        case("s_export_no_stamp", G37, body_plan="no_stamp"),
        case("s_export_bad_row", G37, body_plan="bad_row"),
        case("s_export_other_table", G37, body_plan="other_table"),
        case("s_source_changed_during_export", G37, body_plan="changed_after"),
        # FTP: an existing table.
        case("f_existing_unregistered", FTP, setup=existing(FTP, receipt=False), resource_dir=None),
        case("f_existing_registered", FTP, setup=existing(FTP), resource_dir=None),
        case("f_existing_raw_output_ignored", FTP, setup=existing(FTP), resource_dir=None,
             extra=["--raw-output", "raw/x.tsv"]),
        case("f_existing_different_bytes", FTP, setup=bad_table, resource_dir=None),
        case("f_receipt_without_table", FTP, setup=receipt_only(FTP), resource_dir=None),
        case("f_no_resource_dir", FTP, resource_dir=None),
        case("f_name_conflict", FTP, output="out/" + BASENAME[G37]),
        case("f_resource_in_checkout", FTP, resource_dir="source/res"),
        case("f_resource_is_file", FTP, setup=file_as_folder),
        case("f_lock_held", FTP, hold_lock="out"),
        # FTP: one dump's transfer, stopped at the first failure.
        case("f_part_symlink", FTP, setup=link_raw_part),
        case("f_part_too_long", FTP, setup=long_part),
        case("f_raw_wrong_size", FTP, setup=junk_raw(10)),
        case("f_raw_wrong_sha", FTP, setup=junk_raw(3377160)),
        case("f_range_ignored", FTP, setup=part("CHECKSUMS"), plan={"CHECKSUMS": ["ignore_range"]}),
        case("f_range_refused_shape", FTP, setup=part("CHECKSUMS"), plan={"CHECKSUMS": ["bad_range"]}),
        case("f_wrong_length", FTP, plan={"CHECKSUMS": [{"extra_length": 5, "cut": 100}]}),
        case("f_no_length", FTP, plan={"CHECKSUMS": ["no_length"]}),
        case("f_redirect", FTP, plan={"CHECKSUMS": [{"redirect": "/pub/moved/CHECKSUMS"}]}),
        case("f_404_exhausted", FTP, plan={"CHECKSUMS": [not_found, not_found]}),
        case("f_503_then_ok_then_404", FTP, plan={"CHECKSUMS": [service_failed],
                                                  FTP_FILES[1]: [not_found, not_found]}),
        case("f_cut_body", FTP, plan={"CHECKSUMS": [{"cut": 1_000_000}]}),
        case("f_abrupt_close", FTP, plan={"CHECKSUMS": [{"abrupt": 1_000_000}]}),
        case("f_chunked_reply", FTP, plan={"CHECKSUMS": [{"chunk_cut": 1000}]}),
        case("f_reset_then_resume", FTP, plan={"CHECKSUMS": [{"reset": 1_500_000}],
                                               FTP_FILES[1]: [not_found, not_found]},
             transport=True, volatile=True),
        case("f_complete_part_then_404", FTP, setup=part("CHECKSUMS", size=3377160),
             plan={FTP_FILES[1]: [not_found, not_found]}),
        # FTP: the real projection.
        case("f_full_cold", FTP, attempts=3, timeout=600),
        case("f_full_resume", FTP, attempts=3, timeout=600, setup=resume_seed),
        case("f_full_canonical", FTP, attempts=3, timeout=600, output="source/data/hgnc/" + BASENAME[FTP]),
        case("f_full_retained_dumps", FTP, attempts=3, timeout=600, setup=retained_dumps),
    ]
    return dict(items)


def body_plan(kind, sources, source):
    export = sources.service[source]["export"]
    table = export[:-len(b"[success]\n")]
    lines = table.split(b"\n")
    if kind == "no_stamp":
        return {"export": [{"body": table}]}
    if kind == "bad_row":
        return {"export": [{"body": b"\n".join([lines[0], b"not\ta\tvalid\trow"] + lines[1:])}]}
    if kind == "other_table":
        return {"export": [{"body": b"\n".join([lines[0]] + lines[2:]) + b"[success]\n"}]}
    if kind == "changed_after":
        registry = sources.service[source]["registry"]
        return {"registry": ["ok", {"body": changed_registry(registry)}]}
    raise AssertionError(kind)


def command_line(context, root, spec, program):
    """The fetch arguments, the same for both programs."""
    source = spec["source"]
    base = Path(".") if spec.get("relative") else root
    output = base / spec.get("output", "out/" + BASENAME.get(source, "x.tsv"))
    arguments = ["--source", source, "--output", str(output)]
    raw_default = "raw/" + BASENAME.get(source, "x") + ".raw.tsv"
    raw_output = spec.get("raw_output", None if source == FTP else raw_default)
    if raw_output:
        arguments += ["--raw-output", str(base / raw_output)]
    resource = spec.get("resource_dir", "res" if source == FTP else None)
    if resource:
        arguments += ["--resource-dir", str(base / resource)]
    extra = [str(root / value) if value in ("res", "raw/x.tsv") else value for value in spec.get("extra", [])]
    arguments += extra + ["--attempts", str(spec.get("attempts", 2)), "--timeout", str(spec.get("timeout", 30))]
    if program == "python":
        return [sys.executable, str(root / "source/scripts/biomart_nm_to_enst.py"), "fetch", *arguments]
    return [str(context.binary), "biomart", "fetch", "--helper-root", str(root / "source"), *arguments]


def environment(context, port):
    env = {key: value for key, value in os.environ.items()
           if key.lower() not in ("http_proxy", "https_proxy", "no_proxy", "all_proxy")}
    env.update(https_proxy=f"http://127.0.0.1:{port}", SSL_CERT_FILE=str(context.authority),
               CURL_CA_BUNDLE=str(context.authority), SLURM_JOB_ID="12345", PYTHONDONTWRITEBYTECODE="1")
    return env


def run_measured(command, cwd, env, stdout_path, stderr_path):
    """Run one program; returns exit status, wall seconds and peak memory in KiB."""
    started = time.monotonic()
    with stdout_path.open("wb") as stdout, stderr_path.open("wb") as stderr:
        process = subprocess.Popen(command, cwd=cwd, env=env, stdin=subprocess.DEVNULL, stdout=stdout, stderr=stderr)
        _, status, usage = os.wait4(process.pid, 0)
    process.returncode = os.waitstatus_to_exitcode(status)
    return process.returncode, time.monotonic() - started, usage.ru_maxrss


TIME = re.compile(r"\d{4}-\d\d-\d\dT\d\d:\d\d:\d\d(?:\.\d{1,6})?\+00:00")
TOKEN = re.compile(r"(\.acquisition\.|projection-)[a-z0-9_]{8}")
PRODUCER_FIELD = re.compile(r'("(?:producer_script_path|producer_script_sha256|registration_producer_path|'
                            r'registration_producer_sha256)": )"[^"]*"')


def mask(text, root):
    text = text.replace(str(root), "$ARM")
    text = TIME.sub("$TIME", text)
    text = TOKEN.sub(r"\1$TOKEN", text)
    return PRODUCER_FIELD.sub(r'\1"$PRODUCER"', text)


def inventory(root, volatile):
    """Every file the arm holds, without the copied helper scripts."""
    result = {}
    for path in sorted(root.rglob("*")):
        relative = str(path.relative_to(root))
        if relative.startswith(("source/scripts", "logs")):
            continue
        status = path.lstat()
        item = dict(mode=oct(status.st_mode), links=status.st_nlink)
        if path.is_symlink():
            item["target"] = mask(os.readlink(path), root)
        elif path.is_dir():
            item = dict(mode=oct(status.st_mode))
        elif volatile and path.name.endswith(".part"):
            item.update(bytes="volatile")
        elif path.suffix == ".json" or path.name.endswith(".json.part"):
            item.update(masked_text=hashlib.sha256(mask(path.read_text(), root).encode()).hexdigest())
        else:
            item.update(bytes=status.st_size, sha256=sha256_file(path))
        result[mask(relative, root)] = item
    return result


def failure_of(text, code):
    """The message of a failed run, however the program prints it."""
    text = text.rstrip("\n")
    if code == 0:
        return None
    native = re.fullmatch(r"priva: (\w+): ([\s\S]*)", text)
    if native:
        return native[2]
    if text.startswith("ERROR: "):
        return text[len("ERROR: "):]
    last = text.rsplit("\n", 1)[-1]
    return last.split(": ", 1)[1] if ": " in last else last


def stdout_lines(path, root, transport):
    lines = []
    for line in mask(path.read_text(), root).splitlines():
        if line.startswith("DOWNLOAD "):
            continue
        if line.startswith("Producer "):
            line = "Producer $PRODUCER"
        if transport:
            line = re.sub(r"(failed: ).*(; retrying)", r"\1$TRANSPORT\2", line)
        lines.append(line)
    return lines


def receipts_of(root):
    found = []
    for path in root.rglob("*.json"):
        if not str(path.relative_to(root)).startswith("logs"):
            found.append(path)
    return found


def check_producers(root, program, context, existing):
    """Each new receipt names its own program; the masked text hides these fields."""
    expected_script = str((root / "source/scripts/biomart_nm_to_enst.py").resolve())
    expected_ftp = str(root / "source/scripts/biomart_ftp_nm_to_enst.py")
    problems = []
    for path in receipts_of(root):
        if path in existing:
            continue
        try:
            record = json.loads(path.read_text())
        except ValueError:
            continue
        if not isinstance(record, dict):
            continue
        for path_key, sha_key, script in (("producer_script_path", "producer_script_sha256", expected_script),
                                          ("registration_producer_path", "registration_producer_sha256", expected_ftp)):
            if sha_key not in record:
                continue
            if program == "python":
                wanted_path, wanted_sha = script, sha256_file(script)
            else:
                wanted_path, wanted_sha = NATIVE_PRODUCER, context.builder
            if record[sha_key] != wanted_sha or record.get(path_key, wanted_path) != wanted_path:
                problems.append(f"{path.name}: {record.get(path_key)} {record[sha_key]}")
    return problems


def prepare(root, spec, context):
    root.mkdir(parents=True)
    (root / "logs").mkdir()
    scripts = root / "source/scripts"
    scripts.mkdir(parents=True)
    (root / "source/data/hgnc").mkdir(parents=True)
    for name in ("biomart_nm_to_enst.py", "biomart_ftp_nm_to_enst.py"):
        shutil.copyfile(context.repo / "scripts" / name, scripts / name)
    if "setup" in spec:
        spec["setup"](root, context)
    if "hold_lock" in spec:
        (root / spec["hold_lock"]).mkdir(exist_ok=True)


def run_arm(root, spec, context, program, plan):
    prepare(root, spec, context)
    existing = set(receipts_of(root))
    source = spec["source"]
    logs = root / "logs"
    with Publisher(context.tls, context.sources.content(source), plan) as publisher:
        env = environment(context, publisher.server_address[1])
        lock = None
        if "hold_lock" in spec:
            lock = os.open(root / spec["hold_lock"], os.O_RDONLY)
            fcntl.flock(lock, fcntl.LOCK_EX)
        code, wall, memory = run_measured(command_line(context, root, spec, program), root, env,
                                          logs / "stdout", logs / "stderr")
        if lock is not None:
            os.close(lock)
        requests = list(publisher.requests)
        arm = dict(code=code, wall_seconds=round(wall, 2), peak_rss_kib=memory, requests=requests,
                   stdout=stdout_lines(logs / "stdout", root, spec.get("transport")),
                   failure=mask(failure_of((logs / "stderr").read_text(errors="replace"), code) or "", root),
                   files=inventory(root, spec.get("volatile")),
                   producers=check_producers(root, program, context, existing))
        if code == 0:
            arm.update(after_success(root, spec, context, program, publisher, env))
    return arm


def mapping_path(root, spec):
    source = spec["source"]
    return root / spec.get("output", "out/" + BASENAME[source])


def after_success(root, spec, context, program, publisher, env):
    """A second fetch reuses the table; verify and identity read the first."""
    logs = root / "logs"
    before = len(publisher.requests)
    code, _, _ = run_measured(command_line(context, root, spec, program), root, env,
                              logs / "again.stdout", logs / "again.stderr")
    mapping = mapping_path(root, spec)
    verify = subprocess.run(native_query(context, root, "verify", mapping, spec), capture_output=True, text=True,
                            env=env, cwd=root)
    verify_python = subprocess.run(python_query(root, "verify", mapping, spec), capture_output=True, text=True,
                                   env=env, cwd=root)
    identity = subprocess.run(native_query(context, root, "identity", mapping, spec), capture_output=True, text=True,
                              env=env, cwd=root)
    identity_python = subprocess.run(python_query(root, "identity", mapping, spec), capture_output=True, text=True,
                                     env=env, cwd=root)
    return dict(again_code=code, again_requests=len(publisher.requests) - before,
                again_stdout=stdout_lines(logs / "again.stdout", root, False),
                verify_code=(verify.returncode, verify_python.returncode),
                verify_same=json.loads(verify.stdout or "null") == json.loads(verify_python.stdout or "null"),
                identity_code=(identity.returncode, identity_python.returncode),
                identity_same=identity.stdout == identity_python.stdout,
                table_sha256=sha256_file(mapping))


def native_query(context, root, procedure, mapping, spec):
    command = [str(context.binary), "biomart", procedure, "--helper-root", str(root / "source"), "--mapping", str(mapping)]
    return command + (["--assembly", ASSEMBLY[spec["source"]]] if procedure == "verify" else [])


def python_query(root, procedure, mapping, spec):
    command = [sys.executable, str(root / "source/scripts/biomart_nm_to_enst.py"), procedure, str(mapping)]
    return command + (["--assembly", ASSEMBLY[spec["source"]]] if procedure == "verify" else [])


def comparable_requests(requests, volatile):
    result = []
    for item in requests:
        item = dict(item)
        # A dump request's User-Agent is urllib's own against curl's fixed one.
        if "ensembl_mart_116/" in item["path"]:
            item.pop("user_agent")
        if volatile and item["range"]:
            item["range"] = "present"
        result.append(item)
    return result


def compare(python, native, spec, sources_expected):
    differences = []
    keys = ["code"] if spec.get("usage") else ["code", "stdout", "files"]
    for key in keys:
        if python[key] != native[key]:
            differences.append(key)
    if not spec.get("usage"):
        if comparable_requests(python["requests"], spec.get("volatile")) != \
                comparable_requests(native["requests"], spec.get("volatile")):
            differences.append("requests")
        if not spec.get("transport") and python["failure"] != native["failure"]:
            differences.append("failure")
        if bool(python["failure"]) != bool(native["failure"]):
            differences.append("failure presence")
    if python["producers"] or native["producers"]:
        differences.append("producers")
    if python["code"] == 0 and native["code"] == 0:
        for key in ("again_code", "again_stdout", "again_requests", "table_sha256"):
            if python.get(key) != native.get(key):
                differences.append(key)
        for arm, name in ((python, "python"), (native, "native")):
            if arm["again_code"] != 0 or arm["again_requests"] != 0:
                differences.append(name + " reuse")
            if arm["verify_code"] != (0, 0) or not arm["verify_same"]:
                differences.append(name + " verify")
            if arm["identity_code"] != (0, 0) or not arm["identity_same"]:
                differences.append(name + " identity")
            if arm["table_sha256"] != sources_expected[spec["source"]]:
                differences.append(name + " table sha256")
    return differences


class Context:
    pass


def run_case(name, spec, context):
    spec = dict(spec)
    plan = spec.get("plan")
    if "body_plan" in spec:
        plan = body_plan(spec["body_plan"], context.sources, spec["source"])
    folder = context.work / name
    try:
        python = run_arm(folder / "python", spec, context, "python", plan)
        native = run_arm(folder / "native", spec, context, "native", plan)
    except Exception as failure:
        print(f"{name}: harness failure {failure!r}", flush=True)
        return dict(name=name, differences=[repr(failure)], code=None)
    differences = compare(python, native, spec, context.expected)
    (folder / "comparison.json").write_text(json.dumps(
        dict(name=name, differences=differences, python=python, native=native), indent=1, sort_keys=True))
    print(f"{name}: {'same' if not differences else differences} "
          f"python {python['wall_seconds']}s {python['peak_rss_kib']}KiB "
          f"native {native['wall_seconds']}s {native['peak_rss_kib']}KiB", flush=True)
    for key in ("code", "failure", "stdout", "producers"):
        if python[key] != native[key] or (key == "producers" and python[key]):
            print(f"  {key} python: {str(python[key])[:700]}\n  {key} native: {str(native[key])[:700]}", flush=True)
    if "files" in differences:
        for path in sorted(set(python["files"]) | set(native["files"])):
            if python["files"].get(path) != native["files"].get(path):
                print(f"  file {path}: python {python['files'].get(path)} native {native['files'].get(path)}", flush=True)
    return dict(name=name, differences=differences, code=native["code"],
                python_seconds=python["wall_seconds"], native_seconds=native["wall_seconds"],
                python_rss_kib=python["peak_rss_kib"], native_rss_kib=native["peak_rss_kib"])


def check(args):
    context = Context()
    context.repo = args.repo.resolve()
    context.binary = args.binary.resolve()
    context.raw_dumps = args.raw_dumps
    context.work = args.work_root.resolve()
    context.work.mkdir(parents=True)
    context.tls, context.authority = make_authority(context.work / "authority")
    context.sources = Sources(context.repo, args.raw_dumps)
    context.expected = {source: sha256_file(context.repo / "data/hgnc" / name) for source, name in BASENAME.items()}
    identity = subprocess.run([str(context.binary), "builder-identity", "biomart-fetch"], capture_output=True, text=True,
                              check=True)
    context.builder = json.loads(identity.stdout)["builder_source_sha256"]
    selected = cases()
    if args.only:
        selected = {name: selected[name] for name in args.only}
    if args.skip_full:
        selected = {name: spec for name, spec in selected.items() if name not in FULL}
    # The long cases first, so they overlap with the short ones.
    order = sorted(selected, key=lambda name: name not in FULL)
    with concurrent.futures.ThreadPoolExecutor(args.jobs) as pool:
        results = list(pool.map(lambda name: run_case(name, selected[name], context), order))
    failed = [item for item in results if item["differences"]]
    report = dict(total=len(results), same=len(results) - len(failed), failed=failed,
                  full=[item for item in results if item["name"] in FULL])
    (context.work / "report.json").write_text(json.dumps(dict(report, cases=results), indent=1))
    print(json.dumps(report))
    return 1 if failed else 0


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    commands = parser.add_subparsers(required=True)
    command = commands.add_parser("check")
    command.add_argument("--binary", type=Path, required=True)
    command.add_argument("--repo", type=Path, required=True)
    command.add_argument("--raw-dumps", type=Path, required=True)
    command.add_argument("--work-root", type=Path, required=True)
    command.add_argument("--jobs", type=int, default=4)
    command.add_argument("--only", action="append")
    command.add_argument("--skip-full", action="store_true")
    command.set_defaults(run=check)
    args = parser.parse_args()
    raise SystemExit(args.run(args))


if __name__ == "__main__":
    main()
