#!/usr/bin/env python3
"""Exact, bounded real-cache projection check; not full-cache/VEP acceptance.

Inputs and outputs stay in the authorized local data project. This script contains
no licensed score rows. The independent oracle addresses native columns by name,
rather than using the Rust converter's numeric field offsets.
"""
import argparse
import csv
import hashlib
import json
import os
from pathlib import Path
import shutil
import subprocess
import time

HEADER_SHA256 = "1feb799bed12dae1030e64b03930eab470340848aa4a7e0166072a140e7b0d66"
FIELDS = [f"{effect}_{field}_0" for field in ("delta_score", "dist")
          for effect in ("acceptor_gain", "acceptor_loss", "donor_gain", "donor_loss")]


def sha256(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def main():
    if not __debug__:
        raise SystemExit("Exact validation must run without Python -O/PYTHONOPTIMIZE")
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--prefix-dir", type=Path, required=True)
    parser.add_argument("--context", type=Path, required=True)
    parser.add_argument("--priva-bin", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    args.output_dir.mkdir(exist_ok=False, parents=True)
    with args.context.open() as stream:
        genes = list(csv.DictReader(stream, delimiter="\t"))
    tool_paths = {name: str(Path(shutil.which(name)).resolve()) for name in ("bgzip", "tabix")}
    receipt = {"scope": "512 actual chr1 positive-strand cache-prefix rows only; no inference, whole-cache, indexed VEP or family qualification",
               "converter": {"path": str(args.priva_bin), "sha256": sha256(args.priva_bin)},
               "context": {"path": str(args.context), "sha256": sha256(args.context)},
               "validator_sha256": sha256(__file__), "tools": tool_paths,
               "slurm_job_id": os.environ.get("SLURM_JOB_ID"), "cases": []}
    for kind in ("snv", "indel"):
        source = args.prefix_dir / f"SpliceAI2_GRCh38_{kind}_chr1.prefix256.tsv"
        with source.open() as stream:
            assert hashlib.sha256(stream.readline().encode()).hexdigest() == HEADER_SHA256
            stream.seek(0)
            rows = list(csv.DictReader(stream, delimiter="\t"))
        assert len(rows) == 256
        expected = []
        for row in rows:
            assert row["chrom"] == "chr1" and row["strand"] == "+"
            symbols = sorted({gene["gene_name"] for gene in genes
                              if gene["chrom"] == row["chrom"] and gene["strand"] == row["strand"]
                              and int(gene["gene_start"]) <= int(row["pos"]) <= int(gene["gene_end"])})
            assert symbols, "Native row lacks a verified gene context"
            for symbol in symbols:
                info = "SpliceAI=" + "|".join([row["alt"], symbol] + [row[name] for name in FIELDS])
                expected.append("\t".join([row["chrom"], row["pos"], ".", row["ref"], row["alt"], ".", ".", info]))
        output = args.output_dir / f"{kind}.vcf"
        command = [str(args.priva_bin), "spliceai2-vcf", "--genes", str(args.context),
                   "--chrom", "chr1", "--kind", kind, "--input", str(source), "--output", str(output)]
        started = time.perf_counter()
        result = subprocess.run(command, text=True, capture_output=True, check=True)
        elapsed = time.perf_counter() - started
        (args.output_dir / f"{kind}.stderr.json").write_text(result.stderr)
        counts = json.loads(result.stderr)
        assert counts["input_rows"] == 256 and counts["output_records"] == len(expected)
        content = output.read_text().splitlines()
        actual = [line for line in content if not line.startswith("#")]
        assert actual == expected, f"{kind}: exact full-record/order mismatch"
        assert "##SpliceAI2RecommendedThresholds=0.1,0.25,0.5" in content
        assert "##SpliceAI2Policy=annotation_only_no_ACMG_or_CADD_activation" in content
        compressed = args.output_dir / f"{kind}.vcf.gz"
        with compressed.open("xb") as stream:
            subprocess.run([tool_paths["bgzip"], "-c", str(output)], stdout=stream, check=True)
        subprocess.run([tool_paths["tabix"], "-p", "vcf", str(compressed)], check=True)
        assert subprocess.check_output([tool_paths["tabix"], "-l", str(compressed)], text=True).splitlines() == ["chr1"]
        positions = [int(row["pos"]) for row in rows]
        query = f"chr1:{min(positions)}-{max(positions)}"
        indexed = subprocess.check_output([tool_paths["tabix"], str(compressed), query], text=True).splitlines()
        assert indexed == expected, f"{kind}: exact indexed-record/order mismatch"
        receipt["cases"].append({"kind": kind, "input": str(source), "input_sha256": sha256(source),
                                 "input_rows": len(rows), "output_records": len(actual), "mismatches": 0,
                                 "command": command, "elapsed_seconds": elapsed,
                                 "vcf_sha256": sha256(output), "bgzf_sha256": sha256(compressed),
                                 "index_sha256": sha256(str(compressed) + ".tbi")})
    receipt["status"] = "bounded_actual_projection_passed"
    (args.output_dir / "receipt.json").write_text(json.dumps(receipt, indent=2) + "\n")
    print(json.dumps({"status": receipt["status"], "input_rows": 512, "mismatches": 0}))


if __name__ == "__main__":
    main()
