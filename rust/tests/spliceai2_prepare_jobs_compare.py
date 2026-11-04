#!/usr/bin/env python3
"""Compare `priva spliceai2-prepare prepare --jobs N` with the sequential route.

`--jobs` is a Rust-only option: the shards of each variant kind are projected
and compressed at the same time and the BGZF parts are joined. The compressed
bytes and tabix indexes therefore differ from the sequential pipeline's, so
this script compares what must stay the same: exit status, the decompressed
VCF of each kind (by SHA-256, never printed), and the receipt without the
output byte counts and checksums.

The two runs read a subset of the real source manifest (the chromosomes given
with --chrom, both variant kinds), with vep_gene_annotation supplied from
--vep-annotation as rust/tests/spliceai2_prepare_compare.py does.
"""
import argparse
import gzip
import hashlib
import json
import subprocess
import sys
import tempfile
from pathlib import Path


def sha256_file(path):
    digest = hashlib.sha256()
    with open(path, "rb") as handle:
        for block in iter(lambda: handle.read(1 << 23), b""):
            digest.update(block)
    return digest.hexdigest()


def decompressed_sha256(path):
    digest = hashlib.sha256()
    with gzip.open(path, "rb") as handle:
        for block in iter(lambda: handle.read(1 << 23), b""):
            digest.update(block)
    return digest.hexdigest()


def receipt_view(path, root):
    body = json.loads(Path(path).read_text())
    for key in ("converter", "driver_sha256"):
        body.pop(key, None)
    for output in body.get("outputs", {}).values():
        for key in ("bytes", "sha256", "index_sha256"):
            output.pop(key, None)
    return json.loads(json.dumps(body).replace(str(root), "<ARM>"))


def run(priva, manifest, arm, jobs):
    command = [priva, "spliceai2-prepare", "prepare", "--manifest", str(manifest), "--assembly", "GRCh38",
               "--output-dir", str(arm / "out"), "--receipt", str(arm / "out/adapter.ready.json"),
               "--threads", "1", "--jobs", str(jobs), "--validation-subset"]
    with (arm.parent / f"{arm.name}.log").open("wb") as log:
        status = subprocess.run(command, stdout=log, stderr=subprocess.STDOUT, stdin=subprocess.DEVNULL).returncode
    result = {"exit": status}
    receipt = arm / "out/adapter.ready.json"
    if status == 0 and receipt.exists():
        result["receipt"] = receipt_view(receipt, arm)
        result["vcf"] = {path.name: decompressed_sha256(path) for path in sorted((arm / "out").glob("*.vcf.gz"))}
        result["leftover_parts"] = sorted(path.name for path in (arm / "out").rglob("*.part*"))
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--priva", required=True)
    parser.add_argument("--manifest", required=True, type=Path)
    parser.add_argument("--vep-annotation", required=True, type=Path)
    parser.add_argument("--chrom", action="append", required=True)
    parser.add_argument("--jobs", type=int, default=4)
    parser.add_argument("--work", type=Path)
    options = parser.parse_args()

    data = json.loads(options.manifest.read_text())
    annotation = options.vep_annotation
    data["vep_gene_annotation"] = {"path": str(annotation), "bytes": annotation.stat().st_size,
                                   "sha256": sha256_file(annotation), "release": "GENCODE50"}
    data["files"] = [item for item in data["files"] if item["chrom"] in options.chrom]
    work = Path(tempfile.mkdtemp(prefix="spliceai2_jobs_", dir=options.work))
    subset = work / "subset.json"
    subset.write_text(json.dumps(data, indent=2) + "\n")

    sequential = run(options.priva, subset, work / "sequential", 1)
    parallel = run(options.priva, subset, work / "parallel", options.jobs)
    same = sequential == parallel and sequential["exit"] == 0 and not parallel["leftover_parts"]
    print(json.dumps({"shards": len(data["files"]), "jobs": options.jobs, "same": same,
                      "sequential_exit": sequential["exit"], "parallel_exit": parallel["exit"],
                      "work": str(work)}))
    if not same:
        for key in sorted(set(sequential) | set(parallel)):
            if sequential.get(key) != parallel.get(key):
                print(f"differs: {key}", file=sys.stderr)
    return 0 if same else 1


if __name__ == "__main__":
    sys.exit(main())
