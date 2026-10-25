#!/usr/bin/env python3
"""Verify/reuse one exact GRCh38 BioMart116 FTP projection as a data resource.

This script never selects a PriVA mapping or writes its bundled data/config.
SQL NULL stays NULL in SQLite and serializes as an empty TSV field. Raw row
ordinals determine a new reproducible order; no BioMart service order is claimed.
"""
from __future__ import annotations

import argparse
from collections import Counter
from datetime import datetime, timezone
import fcntl
import gzip
import hashlib
import importlib.util
import json
import os
from pathlib import Path
import re
import sqlite3
import subprocess
import sys
import tempfile
import time
from urllib.request import Request, urlopen
import xml.etree.ElementTree as ET

BASE_URL = "https://ftp.ensembl.org/pub/release-116/mysql/ensembl_mart_116/"
DATASET = "hsapiens_gene_ensembl"
ASSEMBLY = "GRCh38.p14"
HEADERS = ("Gene stable ID", "Transcript stable ID", "HGNC symbol",
           "UniProtKB/Swiss-Prot ID", "RefSeq mRNA ID")
ATTRIBUTES = {
    "ensembl_gene_id": ("gene_id_1020_key", "main", "stable_id_1023", HEADERS[0]),
    "ensembl_transcript_id": ("transcript_id_1064_key", "main", "stable_id_1066", HEADERS[1]),
    "hgnc_symbol": ("gene_id_1020_key", DATASET + "__ox_hgnc__dm", "display_label_1074", HEADERS[2]),
    "uniprotswissprot": ("translation_id_1068_key", DATASET + "__ox_uniprotswissprot__dm", "dbprimary_acc_1074", HEADERS[3]),
    "refseq_mrna": ("transcript_id_1064_key", DATASET + "__ox_refseq_mrna__dm", "dbprimary_acc_1074", HEADERS[4]),
}
# Observed exact publisher CHECKSUMS snapshot, independently retained in raw/.
# Each BSD sum and 1024-byte block count comes from that checksum document.
PINS = {
    "CHECKSUMS": (3377160, None, "da42353dbc969077a61022369c472af4eeea10bcde712599a4c94133d7d8be23"),
    "ensembl_mart_116.sql.gz": (2847327, (42829, 2781), "b615058721edb820e9d4f60da7a21b1bad9088a03fb36047e901b944009d9bf3"),
    "meta_conf__dataset__main.txt.gz": (7126, (12617, 7), "d285b195abdd99db3fe6054d37d5e42f3ed7fe4aa8b80d521c991f12b64f939a"),
    "meta_conf__xml__dm.txt.gz": (34552278, (45435, 33743), "e2b351aea776b65e78a1896a2c4d21d756b8c2771df7db467020f78eb2dcde59"),
    DATASET + "__transcript__main.txt.gz": (35666730, (6815, 34831), None),
    DATASET + "__translation__main.txt.gz": (47541608, (6440, 46428), None),
    DATASET + "__ox_hgnc__dm.txt.gz": (1006750, (4596, 984), None),
    DATASET + "__ox_uniprotswissprot__dm.txt.gz": (1707967, (47451, 1668), None),
    DATASET + "__ox_refseq_mrna__dm.txt.gz": (2357614, (47156, 2303), None),
}
PROJECTIONS = {
    "transcripts": (DATASET + "__transcript__main", ("gene_id_1020_key", "transcript_id_1064_key", "stable_id_1023", "stable_id_1066", "biotype_1064")),
    "main": (DATASET + "__translation__main", ("gene_id_1020_key", "transcript_id_1064_key", "translation_id_1068_key", "stable_id_1023", "stable_id_1066", "biotype_1064")),
    "hgnc": (DATASET + "__ox_hgnc__dm", ("gene_id_1020_key", "display_label_1074")),
    "swissprot": (DATASET + "__ox_uniprotswissprot__dm", ("translation_id_1068_key", "dbprimary_acc_1074")),
    "refseq": (DATASET + "__ox_refseq_mrna__dm", ("transcript_id_1064_key", "dbprimary_acc_1074")),
}
QUERY = """SELECT m.stable_id_1023, m.stable_id_1066, h.display_label_1074,
 u.dbprimary_acc_1074, r.dbprimary_acc_1074
FROM main AS m
LEFT JOIN hgnc AS h ON h.gene_id_1020_key=m.gene_id_1020_key
LEFT JOIN swissprot AS u ON u.translation_id_1068_key=m.translation_id_1068_key
LEFT JOIN refseq AS r ON r.transcript_id_1064_key=m.transcript_id_1064_key
ORDER BY m.__row, h.__row, u.__row, r.__row"""


def require(condition, message):
    if not condition:
        raise ValueError(message)


def sha256(path):
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def atomic_json(path, record):
    part = path.with_name(path.name + ".part")
    require(not path.exists() and not path.is_symlink(), f"Preserve existing receipt: {path}")
    with part.open("x") as handle:
        json.dump(record, handle, indent=2, sort_keys=True)
        handle.write("\n")
    os.replace(part, path)


def verify_raw(path, pin):
    size, bsd, expected_sha = pin
    require(path.is_file() and path.stat().st_size == size, f"Raw size mismatch: {path}")
    digest = sha256(path)
    if expected_sha:
        require(digest == expected_sha, f"Pinned raw SHA256 mismatch: {path}")
    result = {"path": str(path), "bytes": size, "sha256": digest}
    if bsd:
        observed = subprocess.check_output(["sum", "-r", "--", str(path)], text=True).split()
        require(tuple(map(int, observed[:2])) == bsd, f"Publisher BSD checksum mismatch: {path}")
        result["publisher_bsd_sum_blocks"] = list(bsd)
    if path.name.endswith((".gz", ".gz.part")):
        total = 0
        with gzip.open(path, "rb") as handle:
            for block in iter(lambda: handle.read(1024 * 1024), b""):
                total += len(block)
        result["full_gzip_integrity_verified"] = True
        result["uncompressed_bytes"] = total
    return result


def acquire(raw, name, timeout=120):
    path = raw / name
    pin = PINS[name]
    if path.exists() or path.is_symlink():
        record = verify_raw(path, pin)
        print("REUSE", name, record["bytes"], flush=True)
        return {**record, "url": BASE_URL + name, "reused": True}
    part = raw / (name + ".part")
    require(not part.is_symlink(), f"Refuse a symlinked transfer partial: {part}")
    offset = part.stat().st_size if part.exists() else 0
    require(0 <= offset <= pin[0], f"Unexpected partial length: {part}")
    if offset < pin[0]:
        headers = {"Accept-Encoding": "identity"}
        if offset:
            headers["Range"] = f"bytes={offset}-"
        with urlopen(Request(BASE_URL + name, headers=headers), timeout=timeout) as response:
            require(response.url == BASE_URL + name, f"Unexpected download redirect: {response.url}")
            if offset:
                require(response.status == 206 and response.headers.get("Content-Range") == f"bytes {offset}-{pin[0]-1}/{pin[0]}", "Server did not honor exact resume range; partial retained")
            else:
                require(response.status == 200, "Full download must return HTTP200")
            require(int(response.headers["Content-Length"]) == pin[0] - offset, "Publisher transfer length changed")
            with part.open("ab" if part.exists() else "xb") as output:
                received = offset
                reported = time.monotonic()
                for block in iter(lambda: response.read(1024 * 1024), b""):
                    received += len(block)
                    require(received <= pin[0], "Download exceeded pinned size")
                    output.write(block)
                    if time.monotonic() - reported >= 20:
                        print("DOWNLOAD", name, received, "/", pin[0], flush=True)
                        reported = time.monotonic()
    record = verify_raw(part, pin)
    os.replace(part, path)
    record["path"] = str(path)
    print("ACQUIRED", name, record["bytes"], record["sha256"], flush=True)
    return {**record, "url": BASE_URL + name, "reused": False}


def mysql_value(value):
    """Decode MySQL SELECT-INTO-OUTFILE escapes, distinguishing NULL and empty."""
    if value == b"\\N":
        return None
    escapes = {48: 0, 98: 8, 110: 10, 114: 13, 116: 9, 90: 26, 92: 92}
    result = bytearray()
    index = 0
    while index < len(value):
        current = value[index]
        if current == 92:
            index += 1
            require(index < len(value), "Truncated MySQL field escape")
            current = escapes.get(value[index], value[index])
        result.append(current)
        index += 1
    return bytes(result)


def read_schema(path):
    wanted = {table for table, _ in PROJECTIONS.values()}
    wanted.update(("meta_conf__dataset__main", "meta_conf__xml__dm"))
    schemas = {}
    selected = None
    total = 0
    with gzip.open(path, "rt", encoding="latin1") as handle:
        for line in handle:
            total += len(line)
            require(total <= 512_000_000, "Schema exceeds bounded size")
            if line.startswith("CREATE TABLE "):
                name = re.match(r"CREATE TABLE `([^`]+)`", line).group(1)
                selected = name if name in wanted else None
                if selected:
                    require(selected not in schemas, f"Duplicate schema definition: {selected}")
                    schemas[selected] = []
            elif selected and line.startswith("  `"):
                match = re.match(r"  `([^`]+)`\s+(.+)", line)
                require(match is not None, "Cannot parse publisher column schema")
                schemas[selected].append((match.group(1), match.group(2)))
            elif selected and line.startswith(")"):
                require("CHARSET=latin1" in line, f"Unexpected table character set: {selected}")
                selected = None
    require(set(schemas) == wanted, "Required publisher schema tables missing")
    for table, fields in PROJECTIONS.values():
        names = [name for name, _ in schemas[table]]
        require(len(names) == len(set(names)), f"Repeated source column: {table}")
        require(set(fields).issubset(names), f"Required projected fields absent: {table}")
        for name, declaration in schemas[table]:
            if name in fields and name.endswith("_key"):
                require(declaration.startswith("int("), f"Unexpected join key type: {table}.{name}")
    return schemas


def parse_xml(blob, dataset_row):
    require(dataset_row[1:3] == [DATASET, "Human genes (GRCh38.p14)"] and dataset_row[6] == ASSEMBLY, "Publisher dataset/assembly metadata changed")
    root = ET.fromstring(blob)
    require(root.tag == "DatasetConfig" and root.get("dataset") == DATASET and root.get("version") == ASSEMBLY and root.get("datasetID") == dataset_row[0], "Human XML identity differs from dataset metadata")
    require(root.get("displayName") == dataset_row[2] and root.get("interfaces") == "default", "Human dataset display/interface changed")
    records = {}
    for name, (key, table, field, header) in ATTRIBUTES.items():
        matches = [node for node in root.iter("AttributeDescription") if node.get("internalName") == name]
        require(len(matches) == 1, f"Requested attribute is not unique: {name}")
        node = matches[0]
        require(tuple(node.get(k) for k in ("key", "tableConstraint", "field", "displayName")) == (key, table, field, header), f"Requested attribute definition changed: {name}")
        records[name] = dict(node.attrib)
    mains = [node.text for node in root.findall("MainTable")]
    keys = [node.text for node in root.findall("Key")]
    require(mains == [DATASET + "__" + level + "__main" for level in ("gene", "transcript", "translation")], "Unexpected main table hierarchy")
    require(keys == ["gene_id_1020_key", "transcript_id_1064_key", "translation_id_1068_key"], "Unexpected main table key hierarchy")
    return {"dataset_row": dataset_row, "root": dict(root.attrib), "attributes": records, "main_tables": mains, "main_keys": keys, "xml_sha256": hashlib.sha256(blob).hexdigest()}


def metadata(raw, schemas, root):
    dataset_columns = [name for name, _ in schemas["meta_conf__dataset__main"]]
    require(dataset_columns == ["dataset_id_key", "dataset", "display_name", "description", "type", "visible", "version", "modified"], "Dataset metadata column order changed")
    with gzip.open(raw / "meta_conf__dataset__main.txt.gz", "rb") as handle:
        matches = [line.rstrip(b"\n").split(b"\t") for line in handle if line.split(b"\t")[1] == DATASET.encode()]
    require(len(matches) == 1, "Human dataset metadata row absent/nonunique")
    row = [mysql_value(cell).decode("ascii") for cell in matches[0]]
    require([name for name, _ in schemas["meta_conf__xml__dm"]] == ["dataset_id_key", "xml", "compressed_xml", "message_digest"], "XML metadata column order changed")
    found = []
    with gzip.open(raw / "meta_conf__xml__dm.txt.gz", "rb") as handle:
        for line in handle:
            fields = line.rstrip(b"\n").split(b"\t")
            require(len(fields) == 4, "Malformed publisher XML dump row")
            if fields[0] == row[0].encode():
                found.append([mysql_value(cell) for cell in fields])
    require(len(found) == 1, "Human XML row absent/nonunique")
    fields = found[0]
    blob = fields[1]
    require(blob is not None and fields[2] is not None, "Expected raw and compressed human XML copies")
    require(gzip.decompress(fields[2]) == blob, "Raw/compressed human XML copies differ")
    result = parse_xml(blob, row)
    require(result["xml_sha256"] == "a5afb502ad97bada17f9fd4cdc3886be18ba8b237f3b896568793eb0cdd74fa7", "Pinned human XML changed")
    view = root / "hsapiens_gene_ensembl.dataset.xml"
    if view.exists():
        require(view.read_bytes() == blob, "Existing human XML view differs")
    else:
        view.write_bytes(blob)
    return result


def create_projection(db, alias, columns):
    declarations = [f'"{name}" ' + ("INTEGER" if name.endswith("_key") else "TEXT") for name in columns]
    db.execute(f'CREATE TABLE "{alias}" (__row INTEGER PRIMARY KEY, ' + ", ".join(declarations) + ")")


def load_projection(db, raw, alias, table, columns, schema):
    names = [name for name, _ in schema]
    indices = [names.index(name) for name in columns]
    create_projection(db, alias, columns)
    batch, row_number, nulls, blanks = [], 0, Counter(), Counter()
    sql = f'INSERT INTO "{alias}" VALUES (' + ",".join("?" for _ in range(len(columns) + 1)) + ")"
    with gzip.open(raw / (table + ".txt.gz"), "rb") as handle:
        for line in handle:
            row_number += 1
            cells = line.rstrip(b"\n").split(b"\t")
            require(len(cells) == len(names), f"Publisher dump column count differs at {table}:{row_number}")
            projected = []
            for name, index in zip(columns, indices):
                value = mysql_value(cells[index])
                if value is None:
                    nulls[name] += 1
                    projected.append(None)
                elif name.endswith("_key"):
                    require(re.fullmatch(rb"\d+", value) is not None, f"Noninteger key: {table}:{row_number}")
                    projected.append(int(value))
                else:
                    text = value.decode("ascii")
                    require(not any(c in text for c in "\t\n\r\0"), f"Unsafe TSV field: {table}:{row_number}:{name}")
                    if text == "":
                        blanks[name] += 1
                    if name in ("stable_id_1023", "stable_id_1066"):
                        require(re.fullmatch("ENSG[0-9]+" if name == "stable_id_1023" else "ENST[0-9]+", text) is not None, f"Malformed stable ID: {table}:{row_number}")
                    projected.append(text)
            batch.append((row_number, *projected))
            if len(batch) >= 10000:
                db.executemany(sql, batch)
                batch.clear()
    if batch:
        db.executemany(sql, batch)
    require(row_number > 0, f"Empty human table: {table}")
    for name in columns:
        if name.endswith("_key"):
            db.execute(f'CREATE INDEX "{alias}_{name}" ON "{alias}" ("{name}")')
    db.commit()
    print("IMPORTED", alias, row_number, flush=True)
    return {"publisher_table": table, "rows": row_number, "selected_columns": list(columns), "source_columns": names, "sql_null_counts": dict(nulls), "empty_string_counts": dict(blanks)}


def scalar(db, query):
    return db.execute(query).fetchone()[0]


def relationship_audit(db):
    report = {}
    report["transcript_rows"] = scalar(db, "SELECT COUNT(*) FROM transcripts")
    report["main_rows"] = scalar(db, "SELECT COUNT(*) FROM main")
    report["distinct_transcript_keys"] = scalar(db, "SELECT COUNT(DISTINCT transcript_id_1064_key) FROM transcripts")
    require(report["transcript_rows"] == report["distinct_transcript_keys"] and report["transcript_rows"] > 0, "Transcript main key is not unique")
    require(scalar(db, "SELECT COUNT(*) FROM transcripts WHERE transcript_id_1064_key IS NULL OR transcript_id_1064_key<=0 OR gene_id_1020_key IS NULL OR gene_id_1020_key<=0 OR stable_id_1023 IS NULL OR stable_id_1066 IS NULL") == 0, "Transcript main lacks stable IDs or positive keys")
    report["transcripts_missing_from_translation_main"] = scalar(db, "SELECT COUNT(*) FROM transcripts t WHERE NOT EXISTS (SELECT 1 FROM main m WHERE m.transcript_id_1064_key=t.transcript_id_1064_key)")
    require(report["transcripts_missing_from_translation_main"] == 0, "Translation main omits transcripts; no candidate is promoted")
    report["translation_main_missing_or_changed_transcripts"] = scalar(db, "SELECT COUNT(*) FROM main m WHERE NOT EXISTS (SELECT 1 FROM transcripts t WHERE t.transcript_id_1064_key=m.transcript_id_1064_key AND t.gene_id_1020_key=m.gene_id_1020_key AND t.stable_id_1023=m.stable_id_1023 AND t.stable_id_1066=m.stable_id_1066 AND t.biotype_1064 IS m.biotype_1064)")
    require(report["translation_main_missing_or_changed_transcripts"] == 0, "Translation main references absent or inconsistent transcripts")
    for alias, key in (("hgnc", "gene_id_1020_key"), ("swissprot", "translation_id_1068_key"), ("refseq", "transcript_id_1064_key")):
        count = scalar(db, f'SELECT COUNT(*) FROM "{alias}" d WHERE d."{key}">0 AND NOT EXISTS (SELECT 1 FROM main m WHERE m."{key}"=d."{key}")')
        report[alias + "_orphan_positive_key_rows"] = count
        require(count == 0, f"Orphan publisher dimension keys: {alias}")
    report["no_translation_main_rows"] = scalar(db, "SELECT COUNT(*) FROM main WHERE translation_id_1068_key IS NULL OR translation_id_1068_key=0")
    report["noncoding_no_translation_transcripts"] = scalar(db, "SELECT COUNT(DISTINCT transcript_id_1064_key) FROM main WHERE (translation_id_1068_key IS NULL OR translation_id_1068_key=0) AND biotype_1064!='protein_coding'")
    require(report["no_translation_main_rows"] > 0 and report["noncoding_no_translation_transcripts"] > 0, "No untranslated/noncoding transcript witness")
    report["untranslated_with_nonblank_swissprot"] = scalar(db, "SELECT COUNT(*) FROM main m JOIN swissprot u ON u.translation_id_1068_key=m.translation_id_1068_key WHERE (m.translation_id_1068_key IS NULL OR m.translation_id_1068_key=0) AND u.dbprimary_acc_1074 IS NOT NULL AND u.dbprimary_acc_1074!=''")
    require(report["untranslated_with_nonblank_swissprot"] == 0, "Untranslated transcript acquired SwissProt mapping")
    report["expected_output_rows"] = scalar(db, """WITH h AS (SELECT gene_id_1020_key k,COUNT(*) n FROM hgnc GROUP BY 1),
 u AS (SELECT translation_id_1068_key k,COUNT(*) n FROM swissprot GROUP BY 1),
 r AS (SELECT transcript_id_1064_key k,COUNT(*) n FROM refseq GROUP BY 1)
SELECT SUM(COALESCE(h.n,1)*COALESCE(u.n,1)*COALESCE(r.n,1)) FROM main m
LEFT JOIN h ON h.k=m.gene_id_1020_key LEFT JOIN u ON u.k=m.translation_id_1068_key LEFT JOIN r ON r.k=m.transcript_id_1064_key""")
    report["untranslated_output_rows"] = scalar(db, """SELECT COUNT(*) FROM main m
 LEFT JOIN hgnc h ON h.gene_id_1020_key=m.gene_id_1020_key
 LEFT JOIN swissprot u ON u.translation_id_1068_key=m.translation_id_1068_key
 LEFT JOIN refseq r ON r.transcript_id_1064_key=m.transcript_id_1064_key
 WHERE m.translation_id_1068_key IS NULL OR m.translation_id_1068_key=0""")
    require(report["untranslated_output_rows"] >= report["no_translation_main_rows"], "Left joins lost untranslated transcript rows")
    return report


def export_candidate(db, output, audit):
    part = output.with_name(output.name + ".part")
    require(not output.exists() and not output.is_symlink(), "Preserve existing candidate output")
    rows, blanks = 0, Counter()
    with part.open("x", encoding="ascii", newline="\n") as handle:
        handle.write("\t".join(HEADERS) + "\n")
        for values in db.execute(QUERY):
            rows += 1
            for i, value in enumerate(values):
                if value is None or value == "":
                    blanks[HEADERS[i]] += 1
            handle.write("\t".join("" if value is None else value for value in values) + "\n")
    require(rows == audit["expected_output_rows"], "Export row count differs from independent dimension product")
    os.replace(part, output)
    return {"rows": rows, "blank_field_counts": dict(blanks), "bytes": output.stat().st_size, "sha256": sha256(output)}


SOURCE_ID = "ensembl-grch38-mart116-ftp"
BASENAME = "ensembl_GRCh38_mart116_ftp_candidate.tsv"
RELATIVE_PATH = Path("data/hgnc") / BASENAME
EXPECTED_BYTES = 28_785_905
EXPECTED_SHA256 = "dfae01f608b6a9093c2b351ca3fe162af0aacc363fa4a09ba2e4dd5de4850243"
EXPECTED_STATISTICS = {
    "completion_stamp_verified": None, "duplicate_full_rows": 66,
    "duplicate_nm_enst_rows": 90, "genes": 86411, "max_enst_per_nm": 42,
    "nm_accessions": 69679, "nm_rows": 77748, "nm_with_multiple_enst": 4126,
    "non_nm_refseq_accessions": [], "non_nm_refseq_rows": 0, "rows": 691498,
    "rows_without_nm": 613750, "rows_without_refseq": 613750,
    "transcripts": 670670, "unique_nm_enst_pairs": 77658,
}
RAW_SHA256 = {
    "CHECKSUMS": "da42353dbc969077a61022369c472af4eeea10bcde712599a4c94133d7d8be23",
    "ensembl_mart_116.sql.gz": "b615058721edb820e9d4f60da7a21b1bad9088a03fb36047e901b944009d9bf3",
    "meta_conf__dataset__main.txt.gz": "d285b195abdd99db3fe6054d37d5e42f3ed7fe4aa8b80d521c991f12b64f939a",
    "meta_conf__xml__dm.txt.gz": "e2b351aea776b65e78a1896a2c4d21d756b8c2771df7db467020f78eb2dcde59",
    DATASET + "__transcript__main.txt.gz": "7e712dd865bae7a5bff38304b6e56a8bc092074586538d2f7d37c7fdd038a3cd",
    DATASET + "__translation__main.txt.gz": "ea5f780ce729004cff0f701928e2e82adedfea310470e822842da052b4118692",
    DATASET + "__ox_hgnc__dm.txt.gz": "9e5b491ffbbbeaf8e73f107165de31e40d801eb4819b7bf319e60ca58a6ae48f",
    DATASET + "__ox_uniprotswissprot__dm.txt.gz": "0816c3b1ebbb09427f8de7847f0e154bc53e0f9d6159ee626cd2cf5ad35939f9",
    DATASET + "__ox_refseq_mrna__dm.txt.gz": "ac79d31b8fd4a3461b720f7057eba302e6732a5898f90bf070d988eac306d242",
}
PINS = {name: (size, bsd, RAW_SHA256[name]) for name, (size, bsd, _) in PINS.items()}
ORDER_POLICY = "Raw translation-main ordinal, then raw HGNC, SwissProt and RefSeq dimension ordinals; this is a new deterministic order."
NULL_POLICY = "SQL NULL and empty strings remain distinct in SQLite and both serialize as empty TSV fields."


def source():
    """A separate explicitly selectable source; it is never an assembly default."""
    return {"source_id": SOURCE_ID, "basename": BASENAME, "relative_path": RELATIVE_PATH,
            "assembly": ASSEMBLY, "mart_release": 116, "annotation_base_release": 116,
            "aliases": ("hg38", "GRCh38", ASSEMBLY),
            "expected_bytes": EXPECTED_BYTES, "expected_sha256": EXPECTED_SHA256,
            "acquisition_method": "ensembl_mysql_ftp"}


def specification():
    """Immutable method/source facts, independent of registration location/time."""
    return {"source_id": SOURCE_ID, "acquisition_method": "ensembl_mysql_ftp",
            "assembly": ASSEMBLY, "mart_release": 116, "annotation_base_release": 116,
            "source_database": "ensembl_mart_116", "source_root": BASE_URL,
            "dataset": DATASET, "attributes": list(ATTRIBUTES), "header": list(HEADERS),
            "raw_sources": [{"name": name, "url": BASE_URL + name, "bytes": size,
                             "sha256": digest, "publisher_bsd_sum_blocks": list(bsd) if bsd else None}
                            for name, (size, bsd, digest) in PINS.items()],
            "schema_sha256": RAW_SHA256["ensembl_mart_116.sql.gz"],
            "dataset_xml_sha256": "a5afb502ad97bada17f9fd4cdc3886be18ba8b237f3b896568793eb0cdd74fa7",
            "query_sql": QUERY, "query_sha256": hashlib.sha256(QUERY.encode()).hexdigest(),
            "order_policy": ORDER_POLICY, "null_policy": NULL_POLICY,
            "multiplicity_policy": "No DISTINCT, row filtering or deduplication; retain each raw dimension combination.",
            "online_query_executed": False, "original_biomart_export_identity": False,
            "original_biomart_response_order_identity": False,
            "scientific_qualified": False, "production_selected": False}


def table_identity(mapping, validate_table):
    require(mapping.is_file() and mapping.stat().st_size == EXPECTED_BYTES, "Pinned FTP table size differs; existing bytes retained")
    require(sha256(mapping) == EXPECTED_SHA256, "Pinned FTP table SHA256 differs; existing bytes retained")
    observed = validate_table(mapping)
    require(observed == EXPECTED_STATISTICS, "Pinned FTP table schema/order/multiplicity statistics differ")
    require(sha256(mapping) == EXPECTED_SHA256, "FTP mapping changed during verification")
    return observed


def make_provenance(mapping, validate_table, operation):
    observed = table_identity(mapping, validate_table)
    return {**specification(), "bytes": EXPECTED_BYTES, "sha256": EXPECTED_SHA256,
            "statistics": observed, "operation": operation,
            "registered_utc": datetime.now(timezone.utc).isoformat(),
            "registration_producer_path": str(Path(__file__)),
            "registration_producer_sha256": sha256(Path(__file__)),
            "slurm_job_id": os.environ.get("SLURM_JOB_ID"),
            "provenance_basis": "Exact reviewed table snapshot and fixed FTP source/method pins, first extracted by Slurm1515184; registration does not imply a new download or scientific acceptance."}


def verify_record(mapping, record, validate_table, assembly=None):
    if assembly is not None:
        require(assembly in source()["aliases"], f"FTP GRCh38 mapping cannot be used for {assembly}")
    require(isinstance(record, dict), "FTP provenance must be a JSON object")
    for key in ("scientific_qualified", "production_selected", "online_query_executed",
                "original_biomart_export_identity", "original_biomart_response_order_identity"):
        require(record.get(key) is False, "FTP resource provenance must retain explicit unqualified/offline status")
    require(all(record.get(key) == value for key, value in specification().items()), "FTP source/release/assembly/method provenance differs")
    require(not any(key in record for key in ("query_xml", "metadata_before", "metadata_after", "registry_xml", "completion_stamp_verified")),
            "FTP provenance must not claim an online BioMart response or completion stamp")
    require(record.get("bytes") == EXPECTED_BYTES and record.get("sha256") == EXPECTED_SHA256,
            "FTP provenance table identity differs")
    require(record.get("operation") in ("registered_existing_exact_snapshot", "projected_from_verified_ftp_sources"), "Unknown FTP resource operation")
    producer_sha256 = record.get("registration_producer_sha256")
    require(isinstance(producer_sha256, str) and re.fullmatch(r"[0-9a-f]{64}", producer_sha256), "Missing FTP registration producer identity")
    require(record.get("statistics") == table_identity(mapping, validate_table), "FTP provenance statistics differ")
    return record


def verify(mapping, validate_table, assembly=None):
    path = Path(str(mapping) + ".provenance.json")
    with path.open() as handle:
        record = json.load(handle)
    return verify_record(mapping, record, validate_table, assembly)


def install(output, resource_dir, validate_table, timeout=600, attempts=3):
    """Reuse exact known bytes offline, or acquire the minimum pinned FTP dumps."""
    require(1 <= attempts <= 5 and timeout > 0, "Use a positive timeout and 1–5 attempts")
    output = Path(output)
    metadata_path = Path(str(output) + ".provenance.json")
    output.parent.mkdir(parents=True, exist_ok=True)
    fd = os.open(output.parent, os.O_RDONLY)
    try:
        fcntl.flock(fd, fcntl.LOCK_EX | fcntl.LOCK_NB)
        if output.exists() or output.is_symlink():
            if metadata_path.exists() or metadata_path.is_symlink():
                return verify(output, validate_table, ASSEMBLY)
            # The immutable SHA/bytes/statistics pin identifies this already
            # reviewed snapshot. Do not invent a new acquisition or service stamp.
            record = make_provenance(output, validate_table, "registered_existing_exact_snapshot")
            atomic_json(metadata_path, record)
            return verify(output, validate_table, ASSEMBLY)
        require(not metadata_path.exists() and not metadata_path.is_symlink(), "Missing FTP table has existing provenance; restore it explicitly")
        root = Path(resource_dir)
        raw = root / "raw"
        raw.mkdir(parents=True, exist_ok=True)
        for name in PINS:
            for attempt in range(attempts):
                try:
                    acquire(raw, name, timeout)
                    break
                except OSError:
                    if attempt + 1 == attempts:
                        raise
                    time.sleep(min(5 * (attempt + 1), 15))
        checksums = {}
        for line in (raw / "CHECKSUMS").read_text().splitlines():
            cells = line.split()
            if len(cells) == 3 and cells[2] in PINS:
                require(cells[2] not in checksums, "Repeated publisher checksum entry")
                checksums[cells[2]] = tuple(map(int, cells[:2]))
        require(all(checksums.get(name) == bsd for name, (_, bsd, _) in PINS.items() if bsd), "Publisher checksum document differs from fixed FTP pins")
        schemas = read_schema(raw / "ensembl_mart_116.sql.gz")
        metadata(raw, schemas, root)
        work = Path(tempfile.mkdtemp(prefix="projection-", dir=root))
        db = sqlite3.connect(work / "human-projection.sqlite")
        temporary = output.with_name(output.name + ".candidate")
        try:
            db.execute("PRAGMA cache_size=-131072")
            for alias, (table, fields) in PROJECTIONS.items():
                load_projection(db, raw, alias, table, fields, schemas[table])
            audit = relationship_audit(db)
            export_candidate(db, temporary, audit)
            require(db.execute("PRAGMA integrity_check").fetchone() == ("ok",), "FTP projection SQLite integrity failed")
        finally:
            db.close()
        record = make_provenance(temporary, validate_table, "projected_from_verified_ftp_sources")
        # Same-parent hard linking is atomic and refuses replacement. Preserve
        # staged data and source dumps if verification or publication fails.
        os.link(temporary, output)
        temporary.unlink()
        atomic_json(metadata_path, record)
        return verify(output, validate_table, ASSEMBLY)
    finally:
        os.close(fd)

