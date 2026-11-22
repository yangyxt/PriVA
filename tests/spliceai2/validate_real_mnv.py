#!/usr/bin/env python3
"""Bounded local publisher-row regression; no licensed rows are embedded here.

Compare the first 12,000 native indel chr1 rows, including source line 11929,
against an independent column-name oracle. This qualifies projection and indexed
retrieval only; it does not qualify model inference, full VEP or ACMG policy.
"""
import argparse
import csv
import hashlib
import json
import os
from pathlib import Path
import shutil
import subprocess


HEADER_SHA256 = "1feb799bed12dae1030e64b03930eab470340848aa4a7e0166072a140e7b0d66"
FIELDS = [f"{effect}_{field}_0" for field in ("delta_score", "dist")
          for effect in ("acceptor_gain", "acceptor_loss", "donor_gain", "donor_loss")]


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def require(condition, message):
    if not condition:
        raise ValueError(message)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, required=True)
    parser.add_argument("--extraction-receipt", type=Path, required=True)
    parser.add_argument("--context", type=Path, required=True)
    parser.add_argument("--reference-bin", type=Path, required=True)
    parser.add_argument("--priva-bin", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    require(not os.environ.get("PYTHONOPTIMIZE"), "Run without Python optimization")
    data_root = Path("/stornext/snfs170/next-gen/scratch/yangyxt_dt/projects/PriVA_upgrade")
    require(args.output_dir.is_absolute() and args.output_dir.resolve().is_relative_to(data_root / "validation"),
            "Keep licensed validation outputs inside the authorized project validation folder")
    args.output_dir.mkdir(parents=True, exist_ok=False)
    tracked = {str(p.resolve()): digest(p) for p in
               (args.input, args.extraction_receipt, args.context, args.reference_bin, args.priva_bin, Path(__file__))}
    require(digest(args.reference_bin) == "bbc3b7288ed2ea63f1f8455ed87b8a9d08c7bc1fe0a562c1d6134525e2a12dc6",
            "Reference must be the qualified unmodified 225aa4a executable")
    report = {"scope": "12,000 actual indel chr1 rows, exact projection/order and indexed retrieval; no inference/full VEP/scientific acceptance",
              "scientific_qualified": False, "slurm_job_id": os.environ.get("SLURM_JOB_ID"),
              "inputs_sha256": tracked, "passed": False}
    try:
        extraction = json.loads(args.extraction_receipt.read_bytes())
        require(extraction.get("extracted") == str(args.input) and
                extraction.get("extracted_sha256") == digest(args.input) and
                extraction.get("data_rows") == 12000 and extraction.get("offending_source_line") == 11929,
                "Bounded input differs from its byte-preserving extraction receipt")
        manifest_path = Path(extraction["source_manifest"])
        require(digest(manifest_path) == extraction["source_manifest_sha256"] ==
                "fbc6c9fcd1859d6f1605f9343d5659b8420d7397c4152da19daab34c6194e49d",
                "Original publisher source inventory identity differs")
        manifest = json.loads(manifest_path.read_bytes())
        require(extraction["source"] in manifest["files"], "Extracted shard lacks the original publisher identity")
        tracked[str(manifest_path.resolve())] = digest(manifest_path)
        with args.input.open() as stream:
            require(hashlib.sha256(stream.readline().encode()).hexdigest() == HEADER_SHA256,
                    "Native 84-column header identity differs")
            stream.seek(0)
            rows = list(csv.DictReader(stream, delimiter="\t"))
        require(len(rows) == 12000, "Expected the bounded 12,000-row prefix")
        bad = rows[11929 - 2]
        require(tuple(bad[key] for key in ("chrom", "pos", "ref", "alt", "strand")) ==
                ("chr1", "957127", "GG", "AA", "-"), "Original failure locus differs")
        with args.context.open() as stream:
            genes = [gene for gene in csv.DictReader(stream, delimiter="\t") if gene["chrom"] == "chr1"]
        by_strand = {strand: [(int(gene["gene_start"]), int(gene["gene_end"]), gene["gene_name"])
                              for gene in genes if gene["strand"] == strand] for strand in ("+", "-")}
        expected, mnv_expected = [], []
        for row in rows:
            require(row["chrom"] == "chr1", "Wrong input chromosome")
            position = int(row["pos"])
            symbols = sorted({name for start, end, name in by_strand[row["strand"]] if start <= position <= end})
            require(symbols, "Native row lacks a verified scored-gene context")
            for symbol in symbols:
                record = "\t".join([row["chrom"], row["pos"], ".", row["ref"], row["alt"], ".", ".",
                                    "SpliceAI=" + "|".join([row["alt"], symbol] + [row[field] for field in FIELDS])])
                expected.append(record)
                if row is bad:
                    mnv_expected.append(record)
        common = ["spliceai2-vcf", "--genes", str(args.context), "--chrom", "chr1", "--kind", "indel",
                  "--input", str(args.input)]
        original = subprocess.run([str(args.reference_bin), *common, "--output", "-"],
                                  stdout=subprocess.DEVNULL, stderr=subprocess.PIPE, text=True)
        (args.output_dir / "original_failure.log").write_text(original.stderr)
        report["original_exit_code"] = original.returncode
        require(original.returncode == 2 and "native TSV line 11929: indel input requires unequal REF/ALT lengths" in original.stderr,
                "Pinned original did not reproduce the exact recorded failure")
        output = args.output_dir / "indel.vcf"
        candidate = subprocess.run([str(args.priva_bin), *common, "--output", str(output)],
                                   stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
        (args.output_dir / "candidate.stderr.json").write_text(candidate.stderr)
        require(candidate.returncode == 0, "Candidate converter failed")
        counts = json.loads(candidate.stderr)
        require(counts["module"] == "spliceai2-vcf" and counts["chrom"] == "chr1" and counts["kind"] == "indel",
                "Candidate reports the wrong projection module or shard")
        require(counts["input_rows"] == len(rows) and counts["output_records"] == len(expected), "Projection counts differ")
        content = output.read_text().splitlines()
        actual = [line for line in content if not line.startswith("#")]
        if actual != expected:
            mismatches = [i + 1 for i in range(max(len(actual), len(expected)))
                          if i >= len(actual) or i >= len(expected) or actual[i] != expected[i]]
            report["mismatch_output_rows"] = mismatches
            raise ValueError(f"Exact full-record/order mismatches: {len(mismatches)}")
        require("##SpliceAI2RecommendedThresholds=0.1,0.25,0.5" in content, "Threshold metadata changed")
        require("##SpliceAI2Policy=annotation_only_no_ACMG_or_CADD_activation" in content, "Scientific policy changed")
        tools = {name: str(Path(shutil.which(name)).resolve()) for name in ("bgzip", "tabix")}
        compressed = args.output_dir / "indel.vcf.gz"
        with compressed.open("xb") as stream:
            subprocess.run([tools["bgzip"], "-c", str(output)], stdout=stream, check=True)
        subprocess.run([tools["tabix"], "-p", "vcf", str(compressed)], check=True)
        require(subprocess.check_output([tools["tabix"], "-l", str(compressed)], text=True).splitlines() == ["chr1"], "Index contig differs")
        query = subprocess.check_output([tools["tabix"], str(compressed), "chr1:957127-957127"], text=True).splitlines()
        observed = [line for line in query if line.split("\t")[1:5] == ["957127", ".", "GG", "AA"]]
        expected_query = [line for line in expected if line.split("\t")[1:5] == ["957127", ".", "GG", "AA"]]
        require(observed == expected_query and all(line in observed for line in mnv_expected),
                "Indexed retrieval changed scored MNV records, multiplicity or order")
        require(all(digest(Path(p)) == value for p, value in tracked.items()), "Validation input changed")
        report.update(passed=True, input_rows=len(rows), output_records=len(actual), exact_mismatches=0,
                      original_failure_reproduced=True, mnv_records=len(mnv_expected), indexed_mnv_preserved=True,
                      genomic_allele_preserved="chr1:957127 GG>AA; model strand -, no decomposition/reverse-complement",
                      thresholds_preserved=[0.1, 0.25, 0.5], tools=tools,
                      outputs_sha256={p.name: digest(p) for p in (output, compressed, Path(str(compressed) + ".tbi"))})
    except Exception as error:
        report["failure"] = str(error)
        raise
    finally:
        (args.output_dir / "receipt.json").write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps({key: report[key] for key in ("passed", "input_rows", "output_records", "exact_mismatches", "mnv_records", "indexed_mnv_preserved")}))


if __name__ == "__main__":
    main()
