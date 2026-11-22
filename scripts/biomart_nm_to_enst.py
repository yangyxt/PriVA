#!/usr/bin/env python3
"""Acquire or verify explicitly pinned assembly-specific BioMart crosswalks.

GRCh37 uses mart116/annotation75. Explicit VEP116 selects GRCh38 FTP116;
VEP113 and legacy calls without a release retain archive113 on GRCh38.
The legacy constants and function defaults remain GRCh37.
"""
from __future__ import annotations

import argparse
from collections import Counter
from datetime import datetime, timezone
import fcntl
from http.client import IncompleteRead
import hashlib
import importlib.util
import json
import os
from pathlib import Path
import re
import sys
import time
from urllib.parse import urlencode, urlparse
from urllib.request import Request, urlopen
import xml.etree.ElementTree as ET

_FTP_SPEC = importlib.util.spec_from_file_location(
    "priva_biomart_ftp", Path(__file__).with_name("biomart_ftp_nm_to_enst.py"))
FTP = importlib.util.module_from_spec(_FTP_SPEC)
_FTP_SPEC.loader.exec_module(FTP)


SOURCE_ID = "ensembl-grch37-mart116"
ENDPOINT = "https://grch37.ensembl.org/biomart/martservice"
DATASET = "hsapiens_gene_ensembl"
MART = "ENSEMBL_MART_ENSEMBL"
ASSEMBLY = "GRCh37.p13"
BASENAME = "ensembl_GRCh37_mart116_nm_to_enst.tsv"
RELATIVE_PATH = Path("data/hgnc") / BASENAME
# Exact completed snapshot acquired by Slurm job 1514785 on 2026-10-16.
# Keep row order and multiplicity immutable even within the same mart release.
EXPECTED_SNAPSHOT_BYTES = 5_084_679
EXPECTED_SNAPSHOT_SHA256 = "f730640e2411807e7bc38dad141acafc16c3ffca5e4d8ef500976ebd586ddbc7"
GRCH38_SOURCE_ID = "ensembl-grch38-mart113"
GRCH38_ENDPOINT = "https://oct2024.archive.ensembl.org/biomart/martservice"
GRCH38_ASSEMBLY = "GRCh38.p14"
GRCH38_BASENAME = "ensembl_GRCh38_mart113_nm_to_enst.tsv"
GRCH38_RELATIVE_PATH = Path("data/hgnc") / GRCH38_BASENAME
# Exact completed snapshot acquired by Slurm1514907,2026-10-16; duplicates/order stay fixed.
GRCH38_EXPECTED_SNAPSHOT_BYTES = 7_324_483
GRCH38_EXPECTED_SNAPSHOT_SHA256 = "bb5f64c3c5eef06da92c95e50068d716c7e288629733efea6eb34dbab7bc1540"
PINNED_BASENAMES = (BASENAME, GRCH38_BASENAME, FTP.BASENAME)
ATTRIBUTES = ("ensembl_gene_id", "ensembl_transcript_id", "hgnc_symbol",
              "uniprotswissprot", "refseq_mrna")
HEADERS = ("Gene stable ID", "Transcript stable ID", "HGNC symbol",
           "UniProtKB/Swiss-Prot ID", "RefSeq mRNA ID")
QUERY = ('<?xml version="1.0" encoding="UTF-8"?>'
         '<!DOCTYPE Query><Query virtualSchemaName="default" formatter="TSV" '
         'header="1" uniqueRows="0" count="" datasetConfigVersion="0.6" '
         'completionStamp="1"><Dataset name="hsapiens_gene_ensembl" interface="default">'
         + "".join(f'<Attribute name="{name}" />' for name in ATTRIBUTES)
         + '</Dataset></Query>')
PRODUCER_SHA256 = hashlib.sha256(Path(__file__).read_bytes()).hexdigest()


def _source(source_id: str = SOURCE_ID) -> dict:
    # Read snapshot constants at call time to preserve existing callers/fixtures.
    if source_id == FTP.SOURCE_ID:
        return FTP.source()
    if source_id == SOURCE_ID:
        return {"source_id": SOURCE_ID, "endpoint": ENDPOINT, "assembly": ASSEMBLY,
                "basename": BASENAME, "relative_path": RELATIVE_PATH,
                "mart_release": 116, "annotation_base_release": 75,
                "aliases": ("hg19", "GRCh37", ASSEMBLY),
                "expected_bytes": EXPECTED_SNAPSHOT_BYTES, "expected_sha256": EXPECTED_SNAPSHOT_SHA256}
    if source_id == GRCH38_SOURCE_ID:
        return {"source_id": GRCH38_SOURCE_ID, "endpoint": GRCH38_ENDPOINT, "assembly": GRCH38_ASSEMBLY,
                "basename": GRCH38_BASENAME, "relative_path": GRCH38_RELATIVE_PATH,
                "mart_release": 113, "annotation_base_release": 113,
                "aliases": ("hg38", "GRCh38", GRCH38_ASSEMBLY),
                "expected_bytes": GRCH38_EXPECTED_SNAPSHOT_BYTES,
                "expected_sha256": GRCH38_EXPECTED_SNAPSHOT_SHA256}
    raise ValueError(f"Unknown explicitly pinned BioMart source: {source_id}")


def source_for_assembly(assembly: str, vep_release=None) -> dict:
    """Select only a reviewed release; unsupported majors require an explicit path."""
    if vep_release is not None and (isinstance(vep_release, bool) or str(vep_release) not in ("113", "116")):
        raise ValueError(f"No pinned automatic BioMart mapping for VEP {vep_release}; supply an explicit assembly-compatible nm_to_enst_tsv")
    for source_id in (SOURCE_ID, GRCH38_SOURCE_ID):
        source = _source(source_id)
        if assembly in source["aliases"]:
            if source_id == GRCH38_SOURCE_ID and str(vep_release) == "116":
                return _source(FTP.SOURCE_ID)
            return source
    raise ValueError(f"No explicitly pinned BioMart source for assembly: {assembly}")


def source_for_basename(basename: str) -> dict:
    for source_id in (SOURCE_ID, GRCH38_SOURCE_ID, FTP.SOURCE_ID):
        source = _source(source_id)
        if source["basename"] == basename:
            return source
    raise ValueError(f"Unknown pinned BioMart filename: {basename}")


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def provenance_path(mapping: Path) -> Path:
    return Path(str(mapping) + ".provenance.json")


def snapshot_identity(mapping: Path, source_id: str = SOURCE_ID) -> dict:
    source = _source(source_id)
    observed = {"bytes": mapping.stat().st_size, "sha256": sha256(mapping)}
    if observed != {"bytes": source["expected_bytes"], "sha256": source["expected_sha256"]}:
        raise ValueError(f"Pinned BioMart default snapshot size/SHA256 mismatch: {mapping}; retain the source response and recover explicitly")
    return observed


def atomic_json(path: Path, value: dict) -> None:
    temporary = Path(str(path) + ".part")
    if temporary.is_symlink():
        raise ValueError(f"Refusing to replace a symlink used as a staging file: {temporary}")
    with temporary.open("w", encoding="utf-8") as handle:
        json.dump(value, handle, indent=2, sort_keys=True)
        handle.write("\n")
    os.replace(temporary, path)


def _request(parameters: dict, timeout: float, source_id: str = SOURCE_ID):
    source = _source(source_id)
    request = Request(source["endpoint"] + "?" + urlencode(parameters),
                      headers={"User-Agent": "PriVA-pinned-BioMart/1", "Accept-Encoding": "identity"})
    response = urlopen(request, timeout=timeout)
    if urlparse(response.url).hostname != urlparse(source["endpoint"]).hostname:
        response.close()
        raise ValueError(f"BioMart redirected away from the explicitly selected {source['assembly']} host")
    return response


def source_metadata(timeout: float = 60, attempts: int = 3, source_id: str = SOURCE_ID) -> dict:
    """Verify the exact mart and dataset, retaining publisher metadata verbatim."""
    if source_id == FTP.SOURCE_ID:
        raise ValueError("The FTP source uses pinned dump/XML metadata, not an online BioMart registry")
    registry = datasets = None
    for attempt in range(1, attempts + 1):
        try:
            with _request({"type": "registry", "virtualSchema": "default"}, timeout, source_id) as response:
                registry = response.read(1024 * 1024).decode("utf-8")
            with _request({"type": "datasets", "mart": MART}, timeout, source_id) as response:
                datasets = response.read(1024 * 1024).decode("utf-8")
            break
        except (OSError, ValueError) as error:
            if attempt == attempts:
                raise
            print(f"Metadata request {attempt}/{attempts} failed: {error}; retrying the same source", flush=True)
            time.sleep(min(5 * attempt, 15))
    return parse_metadata(registry, datasets, source_id)


def parse_metadata(registry: str, datasets: str, source_id: str = SOURCE_ID) -> dict:
    source = _source(source_id)
    locations = [node.attrib for node in ET.fromstring(registry).iter("MartURLLocation")
                 if node.get("name") == MART]
    release = source["mart_release"]
    required = {"database": f"ensembl_mart_{release}", "displayName": f"Ensembl Genes {release}",
                "host": urlparse(source["endpoint"]).hostname, "path": "/biomart/martservice",
                "serverVirtualSchema": "default"}
    if len(locations) != 1 or any(locations[0].get(key) != value for key, value in required.items()):
        raise ValueError(f"Pinned BioMart source changed; expected {required}, observed {locations}")
    records = [line.split("\t") for line in datasets.splitlines()
               if len(line.split("\t")) >= 5 and line.split("\t")[1] == DATASET]
    if len(records) != 1 or records[0][2] != f"Human genes ({source['assembly']})" or records[0][4] != source["assembly"]:
        raise ValueError(f"Pinned BioMart dataset/assembly changed: {records}")
    return {"registry_record": locations[0], "dataset_record": records[0],
            "registry_xml": registry, "datasets_tsv": datasets}


def validate_table(path: Path, *, raw: bool = False, output: Path | None = None) -> dict:
    """Check every exported row; retain row order and all mapping multiplicity."""
    counts = Counter()
    genes, transcripts, accessions, pairs, full_rows = set(), set(), set(), set(), set()
    ignored_accessions = set()
    mappings = {}
    completed = False
    destination = output.open("wb") if output is not None else None
    try:
        with path.open("rb") as handle:
            header = handle.readline()
            if tuple(header.decode("utf-8").rstrip("\r\n").split("\t")) != HEADERS:
                raise ValueError(f"BioMart table must have exactly the pinned five-column schema: {path}")
            if destination:
                destination.write(header)
            for line_number, line in enumerate(handle, 2):
                if raw and line.rstrip(b"\r\n") == b"[success]":
                    if completed or handle.read().strip():
                        raise ValueError("BioMart completion stamp is not the single terminal record")
                    completed = True
                    break
                cells = line.decode("utf-8").rstrip("\r\n").split("\t")
                if len(cells) != 5 or not re.fullmatch(r"ENSG\d+", cells[0]) or not re.fullmatch(r"ENST\d+", cells[1]):
                    raise ValueError(f"Malformed/truncated BioMart record at line {line_number}: {path}")
                nm = cells[4]
                if nm.startswith("NM_") and not re.fullmatch(r"NM_\d+", nm):
                    raise ValueError(f"BioMart NM accessions must be bare IDs for the unchanged motif reader (line {line_number})")
                if nm and not re.fullmatch(r"[A-Z]{2}_\d+(?:\.\d+)?", nm):
                    raise ValueError(f"Malformed BioMart RefSeq accession at line {line_number}: {path}")
                counts["rows"] += 1
                genes.add(cells[0])
                transcripts.add(cells[1])
                full_rows.add(tuple(cells))
                if nm.startswith("NM_"):
                    counts["nm_rows"] += 1
                    accessions.add(nm)
                    pairs.add((nm, cells[1]))
                    mappings.setdefault(nm, set()).add(cells[1])
                else:
                    counts["rows_without_nm"] += 1
                    if nm:
                        counts["non_nm_refseq_rows"] += 1
                        ignored_accessions.add(nm)
                    else:
                        counts["rows_without_refseq"] += 1
                if destination:
                    destination.write(line)
    finally:
        if destination:
            destination.close()
    if raw and not completed:
        raise ValueError("BioMart export is incomplete: missing terminal [success] completion stamp")
    if not counts["nm_rows"]:
        raise ValueError("BioMart export contains no usable NM-to-ENST mappings")
    return {**{key: counts[key] for key in ("rows", "nm_rows", "rows_without_nm",
                                          "non_nm_refseq_rows", "rows_without_refseq")},
            "non_nm_refseq_accessions": sorted(ignored_accessions),
            "genes": len(genes), "transcripts": len(transcripts),
            "nm_accessions": len(accessions), "unique_nm_enst_pairs": len(pairs),
            "duplicate_full_rows": counts["rows"] - len(full_rows),
            "duplicate_nm_enst_rows": counts["nm_rows"] - len(pairs),
            "nm_with_multiple_enst": sum(len(value) > 1 for value in mappings.values()),
            "max_enst_per_nm": max(map(len, mappings.values())),
            "completion_stamp_verified": completed if raw else None}


def verify_pinned(mapping: Path, assembly: str | None = None) -> dict:
    with provenance_path(mapping).open(encoding="utf-8") as handle:
        record = json.load(handle)
    if not isinstance(record, dict):
        raise ValueError(f"Mapping provenance must be a JSON object: {mapping}")
    source = _source(record.get("source_id"))
    if mapping.name in PINNED_BASENAMES and mapping.name != source["basename"]:
        raise ValueError(f"Known BioMart filename conflicts with its pinned source provenance: {mapping}")
    if source["source_id"] == FTP.SOURCE_ID:
        return FTP.verify_record(mapping, record, validate_table, assembly)
    expected = {"source_id": source["source_id"], "endpoint": source["endpoint"], "mart": MART,
                "mart_release": source["mart_release"], "annotation_base_release": source["annotation_base_release"],
                "assembly": source["assembly"], "dataset": DATASET, "query_xml": QUERY,
                "attributes": list(ATTRIBUTES), "header": list(HEADERS)}
    if any(record.get(key) != value for key, value in expected.items()):
        raise ValueError(f"Mapping provenance does not match the explicitly pinned {source['source_id']} source: {mapping}")
    if assembly is not None and assembly not in source["aliases"]:
        raise ValueError(f"Pinned {source['assembly']} BioMart mapping cannot be selected implicitly for {assembly}; supply an assembly-compatible mapping")
    if not record.get("completion_stamp_verified"):
        raise ValueError(f"Mapping provenance lacks a completed-export receipt: {mapping}")
    before, after = record.get("metadata_before", {}), record.get("metadata_after", {})
    for metadata in (before, after):
        if not isinstance(metadata, dict) or not all(isinstance(metadata.get(key), str)
                                                  for key in ("registry_xml", "datasets_tsv")):
            raise ValueError(f"Mapping provenance metadata must contain publisher XML/TSV strings: {mapping}")
        parsed = parse_metadata(metadata.get("registry_xml", ""), metadata.get("datasets_tsv", ""), source["source_id"])
        if any(metadata.get(key) != parsed[key] for key in ("registry_record", "dataset_record")):
            raise ValueError(f"Mapping provenance metadata is inconsistent: {mapping}")
    if before["registry_record"] != after["registry_record"] or before["dataset_record"] != after["dataset_record"]:
        raise ValueError(f"Mapping source changed during acquisition: {mapping}")
    if record.get("query_sha256") != hashlib.sha256(QUERY.encode()).hexdigest():
        raise ValueError(f"Mapping query identity mismatch: {mapping}")
    actual_identity = snapshot_identity(mapping, source["source_id"])
    if actual_identity != {"bytes": record.get("bytes"), "sha256": record.get("sha256")}:
        raise ValueError(f"Pinned BioMart mapping size/SHA256 mismatch: {mapping}; retain it and recover explicitly")
    observed = validate_table(mapping)
    expected_statistics = record.get("statistics", {})
    if not isinstance(expected_statistics, dict) or set(expected_statistics) != set(observed):
        raise ValueError(f"Mapping provenance is missing full-export statistics: {mapping}")
    if any(observed[key] != value for key, value in expected_statistics.items()
           if key != "completion_stamp_verified"):
        raise ValueError(f"Pinned BioMart mapping statistics mismatch: {mapping}")
    return record


def acquire(output: Path, raw_output: Path | None = None, *, timeout: float = 600, attempts: int = 3,
            source_id: str = SOURCE_ID, resource_dir: Path | None = None) -> dict:
    if source_id == FTP.SOURCE_ID:
        if output.name in PINNED_BASENAMES and output.name != FTP.BASENAME:
            raise ValueError("FTP116 source conflicts with the selected pinned mapping filename")
        # raw_output belongs to the service method. The existing installer may
        # pass it, but an FTP projection never manufactures a service response.
        if (resource_dir is None and not os.path.lexists(output)
                and not os.path.lexists(provenance_path(output))):
            raise ValueError("New FTP acquisition requires --resource-dir outside the checkout")
        if resource_dir is not None and not os.path.lexists(output):
            external_acquisition_path(resource_dir)
        return FTP.install(output, resource_dir,
                           validate_table, timeout, attempts)
    if raw_output is None:
        raise ValueError("Online BioMart acquisition requires its raw_output destination")
    if resource_dir is not None:
        raise ValueError("resource_dir is specific to the FTP source")
    if not os.path.lexists(output):
        external_acquisition_path(raw_output)
    _source(source_id)
    output.parent.mkdir(parents=True, exist_ok=True)
    # Lock the existing directory inode: no persistent lock/state file, and
    # same-directory writers cannot race on the deterministic .part paths.
    descriptor = os.open(output.parent, os.O_RDONLY)
    try:
        try:
            fcntl.flock(descriptor, fcntl.LOCK_EX | fcntl.LOCK_NB)
        except BlockingIOError as error:
            raise ValueError(f"Another mapping writer holds {output.parent}; retry after it finishes") from error
        return _acquire_locked(output, raw_output, timeout=timeout, attempts=attempts, source_id=source_id)
    finally:
        os.close(descriptor)


def _acquire_locked(output: Path, raw_output: Path, *, timeout: float, attempts: int,
                    source_id: str = SOURCE_ID) -> dict:
    source = _source(source_id)
    if output.name in PINNED_BASENAMES and output.name != source["basename"]:
        raise ValueError(f"Selected BioMart source conflicts with the requested pinned output filename: {output}")
    if output.exists() or output.is_symlink():
        record = verify_pinned(output, source["assembly"])
        if record["source_id"] != source_id:
            raise ValueError(f"Existing mapping is not the explicitly requested BioMart source: {output}")
        print(f"Reusing completed pinned mapping: {output} ({record['bytes']} bytes, SHA256 {record['sha256']})", flush=True)
        return record
    metadata_path = provenance_path(output)
    if metadata_path.exists() or metadata_path.is_symlink():
        raise ValueError(f"Mapping is absent but its provenance receipt exists: {metadata_path}; retain it and recover explicitly")
    if raw_output.exists() or raw_output.is_symlink():
        raise ValueError(f"A published raw export already exists: {raw_output}; retain it and select a new output path")
    output.parent.mkdir(parents=True, exist_ok=True)
    raw_output.parent.mkdir(parents=True, exist_ok=True)
    if raw_output == output:
        raise ValueError("The completion-stamped raw response and installed TSV must use different paths")
    print(f"Source {source_id}; Slurm job {os.environ.get('SLURM_JOB_ID', 'none')}; receipt {provenance_path(output)}", flush=True)
    print(f"Producer {Path(__file__).resolve()}; SHA256 {PRODUCER_SHA256}", flush=True)
    before = source_metadata(min(timeout, 60), attempts, source_id)
    atomic_json(Path(str(raw_output) + ".metadata-before.json"),
                {"source_id": source_id, "query_xml": QUERY, "metadata": before,
                 "producer_script_sha256": PRODUCER_SHA256,
                 "producer_script_path": str(Path(__file__).resolve()),
                 "slurm_job_id": os.environ.get("SLURM_JOB_ID"),
                 "retrieval_started_utc": datetime.now(timezone.utc).isoformat()})
    raw_partial = Path(str(raw_output) + ".part")
    table_partial = Path(str(output) + ".part")
    if raw_partial.is_symlink() or table_partial.is_symlink():
        raise ValueError("Refusing to write through a BioMart staging-file symlink")
    for attempt in range(1, attempts + 1):
        try:
            print(f"Requesting {source_id} whole-dataset export, attempt {attempt}/{attempts}", flush=True)
            with _request({"query": QUERY}, timeout, source_id) as response, raw_partial.open("wb") as handle:
                while block := response.read(1024 * 1024):
                    handle.write(block)
            statistics = validate_table(raw_partial, raw=True, output=table_partial)
            break
        except (OSError, IncompleteRead) as error:
            if attempt == attempts:
                raise
            print(f"Export failed: {error}; retrying the same pinned source", flush=True)
            time.sleep(min(10 * attempt, 30))
    after = source_metadata(min(timeout, 60), attempts, source_id)
    if before["registry_record"] != after["registry_record"] or before["dataset_record"] != after["dataset_record"]:
        raise ValueError("BioMart source metadata changed during export; no final mapping installed")
    record = {"source_id": source_id, "endpoint": source["endpoint"], "mart": MART,
              "mart_release": source["mart_release"], "annotation_base_release": source["annotation_base_release"],
              "assembly": source["assembly"],
              "dataset": DATASET, "attributes": list(ATTRIBUTES), "header": list(HEADERS),
              "query_xml": QUERY, "query_sha256": hashlib.sha256(QUERY.encode()).hexdigest(),
              "retrieved_utc": datetime.now(timezone.utc).isoformat(),
              "producer_script_sha256": PRODUCER_SHA256,
              "producer_script_path": str(Path(__file__).resolve()),
              "slurm_job_id": os.environ.get("SLURM_JOB_ID"),
              "metadata_before": before, "metadata_after": after,
              "completion_stamp_verified": True, "statistics": statistics,
              **snapshot_identity(table_partial, source_id),
              "raw_bytes": raw_partial.stat().st_size, "raw_sha256": sha256(raw_partial)}
    os.replace(raw_partial, raw_output)
    os.replace(table_partial, output)
    if is_bundled_mapping(output):
        # The full acquisition record stays beside its external raw export;
        # source data contains portable publisher/method facts only. This is
        # reached only for a newly acquired table, never an existing sidecar.
        runtime_receipt(record, raw_output.parent, output.name)
        record = {key: value for key, value in record.items()
                  if key not in ("retrieved_utc", "producer_script_path", "slurm_job_id")}
        record.update(metadata_kind="priva_bundled_mapping_snapshot", metadata_schema_version=1,
                      acquisition_method="online_biomart",
                      order_policy="Preserve BioMart response row order.",
                      multiplicity_policy="uniqueRows=0; retain duplicate rows and blank RefSeq fields.",
                      scientific_qualified=False, production_selected=False,
                      original_biomart_export_identity=False,
                      original_biomart_response_order_identity=False,
                      online_query_executed=True)
    atomic_json(provenance_path(output), record)
    verify_pinned(output, source["assembly"])
    print(f"Completed export; provenance receipt closed: {provenance_path(output)}", flush=True)
    print(json.dumps({key: record[key] for key in ("source_id", "bytes", "sha256", "statistics")}, sort_keys=True), flush=True)
    return record


def runtime_receipt(record: dict, workspace: Path, basename: str) -> Path:
    """Keep each acquisition's full runtime provenance outside source data."""
    external_acquisition_path(workspace)
    return FTP.runtime_receipt(record, workspace, basename)


def external_acquisition_path(path: Path) -> None:
    """Keep raw exports, FTP dumps and projection databases outside source data."""
    root = Path(__file__).resolve().parents[1]
    if path.resolve().is_relative_to(root.resolve()):
        raise ValueError("BioMart acquisition workspace must be outside the checkout")


def is_bundled_mapping(mapping: Path, source_root: Path | None = None) -> bool:
    """Recognize canonical bundles without trusting symlinked source bytes."""
    root = source_root or Path(__file__).resolve().parents[1]
    logical = Path(os.path.abspath(mapping))
    return (logical.name in PINNED_BASENAMES
            and logical == Path(os.path.abspath(root / "data/hgnc" / logical.name))
            and not logical.is_symlink()
            and logical.resolve() == root.resolve() / "data/hgnc" / logical.name)


def mapping_identity(mapping: Path, inputs: list[Path] | tuple = ()) -> dict:
    metadata = provenance_path(mapping)
    if is_bundled_mapping(mapping):
        if not mapping.is_file() or mapping.stat().st_size == 0:
            raise ValueError(f"Missing or empty bundled mapping: {mapping}")
    elif metadata.exists() or metadata.is_symlink() or mapping.name in PINNED_BASENAMES:
        verify_pinned(mapping)
    return {"mapping_sha256": sha256(mapping), "mapping_bytes": mapping.stat().st_size,
            "provenance_sha256": sha256(metadata) if metadata.exists() else None,
            "anchor_inputs": [{"name": path.name, "bytes": path.stat().st_size,
                               "sha256": sha256(path)} for path in inputs]}


def _artifact_receipt(artifact: Path) -> dict:
    """Verify the existing receipt and artifact without selecting a new mapping."""
    receipt = Path(str(artifact) + ".nm_to_enst.json")
    with receipt.open(encoding="utf-8") as handle:
        recorded = json.load(handle)
    if not isinstance(recorded, dict):
        raise ValueError("output receipt must be a JSON object")
    identity = recorded.get("mapping_identity")
    required = {"mapping_sha256", "mapping_bytes", "provenance_sha256", "anchor_inputs"}
    digest = lambda value: isinstance(value, str) and re.fullmatch(r"[0-9a-f]{64}", value)
    size = lambda value: type(value) is int and value >= 0
    if (not isinstance(identity, dict) or set(identity) != required
            or not digest(identity["mapping_sha256"]) or not size(identity["mapping_bytes"])
            or (identity["provenance_sha256"] is not None and not digest(identity["provenance_sha256"]))
            or not isinstance(identity["anchor_inputs"], list)):
        raise ValueError("output receipt has an invalid mapping identity")
    for anchor in identity["anchor_inputs"]:
        if (not isinstance(anchor, dict) or set(anchor) != {"name", "bytes", "sha256"}
                or not isinstance(anchor["name"], str) or not anchor["name"]
                or not size(anchor["bytes"]) or not digest(anchor["sha256"])):
            raise ValueError("output receipt has an invalid anchor identity")
    if (not artifact.is_file() or not artifact.stat().st_size
            or type(recorded.get("artifact_bytes")) is not int
            or recorded.get("artifact_bytes") != artifact.stat().st_size
            or recorded.get("artifact_sha256") != sha256(artifact)):
        raise ValueError("output content changed or a prior writer left an incomplete artifact")
    return recorded


def check_receipt(artifact: Path, mapping: Path, inputs: list[Path] | tuple = ()) -> None:
    try:
        recorded = _artifact_receipt(artifact)
        if recorded.get("mapping_identity") != mapping_identity(mapping, inputs):
            raise ValueError("mapping content/provenance or HCSeeker/RMC anchors changed")
    except (OSError, ValueError) as error:
        raise ValueError(f"Cannot reuse or overwrite existing {artifact}: missing/stale NM-to-ENST receipt ({error}). Select a new output path or explicitly move the existing artifact and its receipt before rebuilding.") from error


def _managed_artifact_path(*, cache_root: Path, output_dir: Path, relative: Path) -> Path:
    """Resolve an explicit root alias once; reject linked children or escapes."""
    if not cache_root.is_absolute() or not output_dir.is_absolute():
        raise ValueError("Changed mapping outputs require absolute cache_root and output_dir")
    logical_root = Path(os.path.abspath(cache_root))
    root = logical_root.resolve(strict=True)
    if not root.is_dir() or root == root.parent:
        raise ValueError("Changed mapping outputs require a configured cache_root directory")
    logical_output = Path(os.path.abspath(output_dir))
    if logical_output.is_relative_to(logical_root):
        target = root / logical_output.relative_to(logical_root) / relative
    elif logical_output.is_relative_to(root):
        target = logical_output / relative
    else:
        raise ValueError("Changed mapping output_dir escapes cache_root")
    current = root
    for part in target.relative_to(root).parts:
        current = current / part
        if current.is_symlink():
            raise ValueError(f"Linked managed mapping output is not allowed: {current}")
        if current != target and current.exists() and not current.is_dir():
            raise ValueError(f"Managed mapping output parent is not a directory: {current}")
    return target


def changed_mapping_target(*, identity: dict, cache_root: Path, output_dir: Path,
                           assembly: str, kind: str) -> Path:
    """Compute the existing managed namespace without accepting any artifact."""
    aliases = {"hg19": "hg19", "GRCh37": "hg19", "hg38": "hg38", "GRCh38": "hg38"}
    if assembly not in aliases or kind not in ("motifs", "pm1"):
        raise ValueError("Select an explicit supported assembly and motif/PM1 artifact kind")
    assembly = aliases[assembly]
    names = {"motifs": f"AlphaMissense_{assembly}.kde.pkl", "pm1": f"pm1_regions.{assembly}.pkl"}
    key = hashlib.sha256(json.dumps(identity, sort_keys=True, separators=(",", ":")).encode()).hexdigest()
    return _managed_artifact_path(cache_root=cache_root, output_dir=output_dir,
                                  relative=Path("nm-to-enst") / assembly / kind / key / names[kind])


def artifact_route(*, artifact: Path, mapping: Path, cache_root: Path,
                   output_dir: Path, assembly: str, kind: str,
                   inputs: list[Path] | tuple = ()) -> Path:
    """Read-only routing of an intact old output after a validated input change.

    The installer validates the selected mapping/assembly before this call.
    Namespace identity uses the existing receipt fields, without adding VEP
    release, producer revision or mtimes to scientific identity.
    """
    try:
        recorded = _artifact_receipt(artifact)
    except (OSError, ValueError) as error:
        raise ValueError(f"Cannot reuse or overwrite existing {artifact}: missing/stale NM-to-ENST receipt ({error}). Retain the artifact and recover its receipt explicitly.") from error
    if not mapping.is_file() or not mapping.stat().st_size:
        raise ValueError(f"Missing or empty selected mapping: {mapping}")
    identity = mapping_identity(mapping, inputs)
    if recorded["mapping_identity"] == identity:
        return artifact
    target = changed_mapping_target(identity=identity, cache_root=cache_root,
                                    output_dir=output_dir, assembly=assembly, kind=kind)
    receipt = Path(str(target) + ".nm_to_enst.json")
    if receipt.is_symlink():
        raise ValueError(f"Linked managed mapping receipt is not allowed: {receipt}")
    if target.exists() or receipt.exists():
        if target.is_file() and target.stat().st_nlink != 1:
            raise ValueError(f"Shared managed mapping artifact is not allowed: {target}")
        check_receipt(target, mapping, inputs)
    return target


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    commands = parser.add_subparsers(dest="command", required=True)
    fetch = commands.add_parser("fetch", help="Explicit same-source acquisition; completed pinned output is reused offline")
    fetch.add_argument("--source", choices=(SOURCE_ID, GRCH38_SOURCE_ID, FTP.SOURCE_ID), required=True)
    fetch.add_argument("--output", type=Path, required=True)
    fetch.add_argument("--raw-output", type=Path, help="Raw response destination for online BioMart sources")
    fetch.add_argument("--resource-dir", type=Path, help="External FTP raw dumps/projection workspace; required for new acquisition")
    fetch.add_argument("--timeout", type=float, default=600)
    fetch.add_argument("--attempts", type=int, default=3)
    verify = commands.add_parser("verify")
    verify.add_argument("mapping", type=Path)
    verify.add_argument("--assembly")
    identity = commands.add_parser("identity")
    identity.add_argument("mapping", type=Path)
    identity.add_argument("--input", type=Path, action="append", default=[])
    for name in ("receipt-check", "receipt-write"):
        receipt = commands.add_parser(name)
        receipt.add_argument("artifact", type=Path)
        receipt.add_argument("mapping", type=Path)
        receipt.add_argument("--input", type=Path, action="append", default=[])
        if name == "receipt-write":
            receipt.add_argument("--expected-identity", required=True)
    route = commands.add_parser("artifact-route", help="Read-only managed routing after an intact output's mapping identity changes")
    route.add_argument("--artifact", type=Path, required=True)
    route.add_argument("--mapping", type=Path, required=True)
    route.add_argument("--cache-root", type=Path, required=True)
    route.add_argument("--output-dir", type=Path, required=True)
    route.add_argument("--assembly", choices=("hg19", "GRCh37", "hg38", "GRCh38"), required=True)
    route.add_argument("--kind", choices=("motifs", "pm1"), required=True)
    route.add_argument("--input", type=Path, action="append", default=[])
    route.add_argument("--print-path", action="store_true")
    args = parser.parse_args()
    if args.command == "fetch" and args.source != FTP.SOURCE_ID and args.raw_output is None:
        parser.error("--raw-output is required for online BioMart sources")
    try:
        if args.command == "fetch":
            if not 1 <= args.attempts <= 5 or args.timeout <= 0:
                raise ValueError("Use 1-5 bounded attempts and a positive timeout")
            acquire(args.output, args.raw_output, timeout=args.timeout, attempts=args.attempts,
                    source_id=args.source, resource_dir=args.resource_dir)
        elif args.command == "verify":
            print(json.dumps(verify_pinned(args.mapping, args.assembly), sort_keys=True))
        elif args.command == "identity":
            print(json.dumps(mapping_identity(args.mapping, args.input), sort_keys=True))
        elif args.command == "receipt-check":
            check_receipt(args.artifact, args.mapping, args.input)
        elif args.command == "artifact-route":
            target = artifact_route(artifact=args.artifact, mapping=args.mapping,
                                    cache_root=args.cache_root, output_dir=args.output_dir,
                                    assembly=args.assembly, kind=args.kind, inputs=args.input)
            print(str(target) if args.print_path else json.dumps({"artifact": str(target)}, sort_keys=True))
        else:
            if not args.artifact.is_file() or not args.artifact.stat().st_size:
                raise ValueError(f"Scientific command did not create a nonempty output; no mapping receipt written: {args.artifact}")
            expected = json.loads(args.expected_identity)
            if expected != mapping_identity(args.mapping, args.input):
                raise ValueError("Mapping content/provenance changed during scientific generation; output has no valid mapping receipt")
            atomic_json(Path(str(args.artifact) + ".nm_to_enst.json"),
                        {"mapping_identity": expected, "artifact_bytes": args.artifact.stat().st_size,
                         "artifact_sha256": sha256(args.artifact)})
    except (OSError, ValueError, ET.ParseError, IncompleteRead) as error:
        print(f"ERROR: {error}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
