#!/usr/bin/env python3
"""Run PriVA filtering and ACMG assignment on bundled, prepared demo fixtures.

This is deliberately separate from raw annotation: it never runs VEP, CADD
models, downloads, or the installer, and never grants scientific qualification.
Use run --assembly hg19|hg38 [--output-dir DIR], prepare, or validate --config YAML.
Relative resource paths in a smoke configuration are relative to the checkout.
"""
from __future__ import annotations

import argparse
from collections import Counter
import csv
import gzip
import math
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile

import yaml

ROOT = Path(__file__).resolve().parents[1]
ASSEMBLIES = {"hg19": "hg19", "GRCh37": "hg19", "hg38": "hg38", "GRCh38": "hg38"}
REQUIRED_RESOURCES = (
    "input_vcf", "hub_cadd_file", "clinvar_aa_stat", "clinvar_splice_stat",
    "clinvar_patho_af_stat", "clinvar_patho_exon_af_stat", "clinvar_gene_stat",
    "alphamissense_tranx_domain_map", "alphamissense_intolerant_domains",
    "all_intolerant_domains", "pm1_regions_pkl", "interpro_mapping_pickle",
    "clingen_map", "gene_dosage_sensitivity", "hpo_assertions",
    "hpo_condition_mechanism_json", "gene_nonlof_mechanism_json",
)


def resource(value: str) -> Path:
    path = Path(value).expanduser()
    return path.resolve() if path.is_absolute() else (ROOT / path).resolve()


def require_file(path: Path) -> None:
    if not path.is_file() or path.stat().st_size == 0:
        raise ValueError(f"Missing or empty smoke resource: {path}")


def load_config(path: Path) -> dict:
    with path.open() as stream:
        config = yaml.safe_load(stream)
    if not isinstance(config, dict) or config.get("smoke_fixture_only") is not True:
        raise ValueError("Use a smoke configuration with smoke_fixture_only: true.")
    if config.get("ped_file"):
        raise ValueError("The bundled smoke cases are sites-only; ped_file must be empty.")
    return config


def vcf_rows(path: Path) -> list[tuple[str, str, str, str]]:
    rows = []
    fields = None
    with gzip.open(path, "rt") as stream:
        for line in stream:
            if line.startswith("##INFO=<ID=CSQ"):
                fields = line.split("Format: ", 1)[1].split('"', 1)[0].split("|")
            elif not line.startswith("#"):
                values = line.rstrip("\n").split("\t")
                rows.append((values[0], values[1], values[3], values[4]))
                csq = next((x[4:] for x in values[7].split(";") if x.startswith("CSQ=")), "")
                if fields is None or not csq:
                    raise ValueError("Smoke input must contain completed CSQ annotations.")
                if any(len(row.split("|")) != len(fields) for row in csq.split(",")):
                    raise ValueError("Smoke input has a CSQ row/header width mismatch.")
    if not rows:
        raise ValueError("Smoke VCF contains no records.")
    return rows


def input_check(config: dict) -> None:
    for key in REQUIRED_RESOURCES:
        if not config.get(key):
            raise ValueError(f"Missing smoke configuration key: {key}")
        require_file(resource(config[key]))
    variant_path = resource(config["input_vcf"])
    require_file(Path(str(variant_path) + ".tbi"))
    require_file(ROOT / "data/repeats" / config["repeat_region_file_name"])
    variants = set(vcf_rows(variant_path))
    with resource(config["hub_cadd_file"]).open() as stream:
        reader = csv.DictReader(stream, delimiter="\t")
        needed = {"#Chrom", "Pos", "Ref", "Alt", "FeatureID", "PHRED"}
        if not needed.issubset(reader.fieldnames or ()):
            raise ValueError("Prepared CADD fixture needs allele, FeatureID, and PHRED columns.")
        covered = {(r["#Chrom"], r["Pos"], r["Ref"], r["Alt"]) for r in reader}
    # Compare the existing chromosome spellings; never silently rewrite alleles.
    if not variants.issubset(covered):
        raise ValueError("Prepared CADD fixture does not cover all input alleles.")


def prepare(args) -> Path:
    assembly = ASSEMBLIES[args.assembly]
    source = args.config or ROOT / f"example/configs/smoke_demo20_{assembly}.yaml"
    config = load_config(Path(source).resolve())
    if ASSEMBLIES.get(config.get("assembly")) != assembly:
        raise ValueError("Requested assembly and smoke configuration disagree.")
    for key in REQUIRED_RESOURCES:
        if config.get(key):
            config[key] = str(resource(config[key]))
    input_check(config)
    if args.output_dir:
        folder = args.output_dir.expanduser().resolve()
        if folder == ROOT or ROOT in folder.parents or folder in ROOT.parents:
            raise ValueError("Smoke output must be outside the source checkout.")
        if folder.exists() and any(folder.iterdir()):
            raise ValueError(f"Smoke output directory is not empty; retained: {folder}")
        folder.mkdir(parents=True, exist_ok=True)
    else:
        folder = Path(tempfile.mkdtemp(prefix=f"priva-smoke-{assembly}-")).resolve()
    for name in ("input", "output", "tmp"):
        (folder / name).mkdir()
    staged = folder / "input" / f"demo20.{assembly}.vcf.gz"
    shutil.copy2(config["input_vcf"], staged)
    shutil.copy2(config["input_vcf"] + ".tbi", str(staged) + ".tbi")
    config.update(input_vcf=str(staged), control_vcf=str(staged),
                  output_dir=str(folder / "output"), tmp_dir=str(folder / "tmp"),
                  threads=args.cores, threads_per_fam=min(args.cores, int(config.get("threads_per_fam", 1))))
    path = folder / "smoke.config.yaml"
    with path.open("x") as stream:
        yaml.safe_dump(config, stream, sort_keys=False)
    return path


def validate(args) -> None:
    import pandas as pd

    config = load_config(args.config)
    base = Path(config["input_vcf"]).name.removesuffix(".vcf.gz")
    prefix = Path(config["output_dir"]) / f"{base}.anno.filtered"
    table = Path(str(prefix) + ".tsv")
    matrix = Path(str(prefix) + ".acmg.tsv")
    prepared = Path(config["output_dir"]) / f"{base}.anno.prepared.tsv"
    for path in (table, matrix, prepared, Path(str(prefix) + ".vcf.gz"), Path(str(prefix) + ".cadd.tsv")):
        require_file(path)
    with table.open() as stream:
        reader = csv.DictReader(stream, delimiter="\t")
        required = {"ACMG_criteria", "ACMG_class", "ACMG_quant_score", "sort_index", "variant_rank"}
        if not required.issubset(reader.fieldnames or ()):
            raise ValueError("Ranked smoke table is missing its ACMG/ranking columns.")
        rows = list(reader)
    if not rows:
        raise ValueError("Smoke prioritization produced no transcript rows.")
    for row in rows:
        if any(not math.isfinite(float(row[key])) for key in ("ACMG_quant_score", "sort_index", "variant_rank")):
            raise ValueError("Nonfinite score/rank in smoke result.")
    with matrix.open() as stream:
        matrix_rows = list(csv.DictReader(stream, delimiter="\t"))
    # ACMG emits its matrix before the existing protein-coding ranking filter.
    # Check both populations against the actual prepared input, with its source
    # row identity/multiplicity, rather than assuming the two output sizes match.
    before = pd.read_table(prepared, low_memory=False).drop_duplicates()
    protein = before.loc[before["BIOTYPE"].eq("protein_coding")]
    if len(matrix_rows) != len(before) or len(rows) != len(protein):
        raise ValueError("Smoke matrix/ranking row counts do not match the prepared transcript populations.")
    keys = ["chrom", "pos", "ref", "alt", "Feature"]
    expected = Counter(protein[keys].astype(str).itertuples(index=False, name=None))
    observed = Counter(tuple(row[key] for key in keys) for row in rows)
    if expected != observed:
        raise ValueError("Ranked smoke allele/FeatureID rows or their multiplicities changed.")
    count = len(vcf_rows(Path(str(prefix) + ".vcf.gz")))
    print(f"PASS {config['assembly']}: {count} retained variants; {len(rows)} ranked and {len(matrix_rows)} ACMG rows. "
          f"Prepared fixtures only; VEP and CADD model calculation not exercised. Outputs: {config['output_dir']}")


def run(args) -> None:
    for tool in ("bash", "python", "yq", "bcftools", "tabix", "snakemake"):
        if not shutil.which(tool):
            raise ValueError(f"Activate the PriVA environment; missing executable: {tool}")
    path = prepare(args)
    print(f"Prepared-fixture smoke test: {path}", flush=True)
    subprocess.run(["snakemake", "--snakefile", str(ROOT / "example/smoke.Snakefile"),
                    "--configfile", str(path), "--directory", str(path.parent),
                    "--cores", str(args.cores), "--rerun-incomplete"], check=True)
    validate(argparse.Namespace(config=path))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    commands = parser.add_subparsers(dest="command", required=True)
    for name, function in (("run", run), ("prepare", prepare)):
        command = commands.add_parser(name)
        command.add_argument("--assembly", choices=tuple(ASSEMBLIES), required=True)
        command.add_argument("--config", type=Path)
        command.add_argument("--output-dir", type=Path)
        command.add_argument("--cores", type=int, default=2)
        command.set_defaults(function=function)
    command = commands.add_parser("validate")
    command.add_argument("--config", type=Path, required=True)
    command.set_defaults(function=validate)
    args = parser.parse_args()
    if getattr(args, "cores", 1) < 1:
        parser.error("--cores must be positive")
    try:
        result = args.function(args)
        if result:
            print(result)
    except (ValueError, OSError, subprocess.CalledProcessError) as error:
        parser.exit(1, f"Smoke test failed: {error}\n")


if __name__ == "__main__":
    main()
