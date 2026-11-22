#!/usr/bin/env python3
"""Fail closed if the reference software environment is incomplete or drifted."""
import argparse
import hashlib
import importlib
import importlib.metadata
import json
import os
from pathlib import Path
import re
import shutil
import subprocess
import sys


def check(spec):
    import yaml
    config = yaml.safe_load(spec.read_text())
    prefix = Path(sys.prefix).resolve()
    records = {}
    for item in (prefix / "conda-meta").glob("*.json"):
        package = json.loads(item.read_text())
        records[package["name"]] = package
    failures = []
    for dependency in config["dependencies"]:
        match = re.fullmatch(r"([\w-]+)(?:=([\w.]+))?", dependency)
        if not match:  # Solver enforces comparison operators, e.g. >=.
            continue
        name, expected = match.groups()
        found = records.get(name)
        if not found:
            failures.append(f"Missing Conda package: {name}")
        elif expected and found["version"] != expected:
            failures.append(f"{name}: expected {expected}, found {found['version']}")
    modules = (
        "numpy", "pandas", "scipy", "scipy.stats", "sklearn", "statsmodels",
        "matplotlib", "polars", "pyarrow", "pysam", "pybedtools", "pyfaidx",
        "diptest", "requests", "yaml", "jsonschema", "openpyxl", "gffutils",
    )
    for name in modules:
        try:
            importlib.import_module(name)
        except Exception as error:
            failures.append(f"Cannot import {name}: {error}")
    executables = {}
    for name in ("vep", "vep_install", "python", "snakemake", "bcftools", "samtools",
                 "bedtools", "bgzip", "tabix", "yq", "parallel", "sambamba", "seqtk", "tsv-pretty"):
        executable = shutil.which(name)
        if not executable or not Path(executable).is_relative_to(prefix):
            failures.append(f"Missing environment-local executable: {name}")
        else:
            executables[name] = executable
    for command in (["vep", "--help"], ["samtools", "--version"], ["bcftools", "--version"],
                    ["snakemake", "--version"], ["perl", "-MBio::DB::HTS", "-e", "1"]):
        try:
            run = subprocess.run(command, text=True, capture_output=True, timeout=60)
            if run.returncode:
                failures.append(f"{' '.join(command)}: {run.stderr[:2500]}")
        except (OSError, subprocess.TimeoutExpired) as error:
            failures.append(f"{' '.join(command)}: {error}")
    if failures:
        raise RuntimeError("\n".join(failures))
    return {
        "schema_version": 1, "status": "software_checked",
        "scope": "Software only; annotation resources and family parity are not established",
        "prefix": str(prefix), "python": sys.version,
        "spec_sha256": hashlib.sha256(spec.read_bytes()).hexdigest(),
        "packages": {name: record["version"] for name, record in sorted(records.items())},
        "executables": executables,
        "slurm_job_id": os.environ.get("SLURM_JOB_ID"),
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--spec", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    result = check(args.spec)
    temporary = args.output.with_name(args.output.name + ".partial")
    temporary.write_text(json.dumps(result, indent=2) + "\n")
    temporary.replace(args.output)
    print(f"Validated PriVA software: {result['prefix']}")


if __name__ == "__main__":
    main()
