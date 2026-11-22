#!/usr/bin/env python3
"""Runs the family steps both ways on a test family and compares them.

Old: the shell functions the analysis tasks ran before phase 4, as they ran
them: `scripts/filtration_vcf_per_fam.sh` for filtration,
`scripts/prioritization_vcf_per_fam.sh` for the whole prioritization, and the
task's step wrapper (tests/fixtures/family-steps-952811c/prioritization_step.sh)
for one prioritization step. Phase 7 removed the scripts, so the old route
takes them from a checkout of commit 841e200 (--old-root, for example
`git worktree add OLD 841e200`). New: `priva task step`, the Rust steps
`priva task run` calls now.

The test family is data/test_example/TEST_FAM.filtered.vcf.gz with a made-up
AF_grpmax_joint for each record and a CADD table made from its records (with
rows of other variants). Filtration runs the real `priva pedigree-filter`.
Every other priva command is a stand-in that logs its arguments and writes
the file the command writes, so a run compares the commands, not the
science, and the ACMG step needs no installed resources (they are empty
files here).

For each case both routes start from the same files in their own folder and
must give the same exit status, the same priva commands (folder names
replaced by ROUTE), the same files in result/, the same records and header of
each VCF (the ##bcftools_ lines, which hold temporary names and dates, left
out) and the same bytes of every other file but an index. Exit 1 on any
difference.

  python3 tests/family_step_compare.py --priva rust/target/release/priva \\
      --old-root OLD --work W [--case NAME]...
"""
import argparse
import gzip
import json
import os
import random
import shutil
import stat
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
FIXTURE = ROOT / "tests" / "fixtures" / "family-steps-952811c" / "prioritization_step.sh"
FAMILY_VCF = ROOT / "data" / "test_example" / "TEST_FAM.filtered.vcf.gz"
FAMILY_PED = ROOT / "data" / "test_example" / "TEST_FAM.ped"

FAKE_PRIVA = r'''#!/usr/bin/env python3
import json, os, sys
args = sys.argv[1:]
if args == ["--version"]:
    print("priva 0.0.0-stand-in")
    sys.exit(0)
with open(os.environ["FAKE_PRIVA_LOG"], "a") as log:
    log.write(json.dumps(args) + "\n")
if args and args[0] == "pedigree-filter":
    os.execv(os.environ["REAL_PRIVA"], [os.environ["REAL_PRIVA"]] + args)
def value(name):
    return args[args.index(name) + 1] if name in args else None
if args and args[0] == "combine-annotations":
    with open(value("-o"), "w") as out:
        out.write("Chr\tStart\tRef\tAlt\tFeature\n22\t100\tA\tG\tENST1\n")
elif args and args[0] in ("splicing-analysis", "utr-anno-interpret"):
    table = value("--anno_table") or value("--variants_table")
    with open(table, "a") as out:
        out.write("# " + args[0] + "\n")
elif args and args[0] == "acmg-assign":
    table = value("--anno_table")
    with open(table[:-4] + ".acmg.tsv", "w") as out:
        out.write("variant_id\tPVS1\n22:100:A:G\t1\n")
'''

RESOURCES = [
    "clinvar_aa_stat", "clinvar_splice_stat", "clinvar_pm5_index",
    "all_intolerant_domains", "alphamissense_intolerant_domains",
    "alphamissense_tranx_domain_map", "pm1_regions_pkl", "clinvar_patho_af_stat",
    "clinvar_patho_exon_af_stat", "interpro_mapping_pickle", "clingen_map",
    "gene_dosage_sensitivity", "hpo_condition_mechanism_json",
    "gene_nonlof_mechanism_json", "relevant_gene_list", "dispensable_gene_list",
    "cds_fasta_file", "hpo_assertions",
]


def run(command, **kwargs):
    return subprocess.run(command, check=True, **kwargs)


def build_family(work: Path) -> dict:
    """The test family's VCF (with AF_grpmax_joint), CADD table and PED."""
    family = work / "family"
    if (family / "done").exists():
        return json.loads((family / "done").read_text())
    shutil.rmtree(family, ignore_errors=True)
    family.mkdir(parents=True)
    rng = random.Random(4)
    records = run(["bcftools", "query", "-f", "%CHROM\t%POS\t%REF\t%ALT\n", str(FAMILY_VCF)],
                  capture_output=True, text=True).stdout.splitlines()
    annotations = family / "af.tsv"
    with open(annotations, "w") as out:
        for line in records:
            chrom, pos, ref, alt = line.split("\t")
            roll = rng.random()
            if roll < 0.15:
                continue  # no AF_grpmax_joint
            af = 0.5 * rng.random() ** 3
            out.write(f"{chrom}\t{pos}\t{ref}\t{alt}\t{af:.5f}\n")
    run(["bgzip", "-f", str(annotations)])
    run(["tabix", "-f", "-s1", "-b2", "-e2", str(annotations) + ".gz"])
    header = family / "af.hdr"
    header.write_text('##INFO=<ID=AF_grpmax_joint,Number=1,Type=Float,Description="test">\n')
    vcf = family / "TEST.anno.vcf.gz"
    run(["bcftools", "annotate", "-a", str(annotations) + ".gz", "-h", str(header),
         "-c", "CHROM,POS,REF,ALT,AF_grpmax_joint", "-Oz", "-o", str(vcf), str(FAMILY_VCF)])
    run(["tabix", "-f", "-p", "vcf", str(vcf)])
    cadd = family / "TEST.anno.cadd.tsv"
    with open(cadd, "w") as out:
        out.write("## CADD GRCh38-v1.7 (c) University of Washington\n")
        out.write("#Chrom\tPos\tRef\tAlt\tRawScore\tPHRED\n")
        for line in records:
            chrom, pos, ref, alt = line.split("\t")
            if rng.random() < 0.9:
                out.write(f"{chrom.removeprefix('chr')}\t{pos}\t{ref}\t{alt}\t"
                          f"{rng.random():.4f}\t{rng.random() * 40:.3f}\n")
            if rng.random() < 0.05:
                out.write(f"{chrom.removeprefix('chr')}\t{pos}\t{ref}\tN\t0.1\t1.0\n")
    empty = family / "TEST.empty.vcf.gz"
    run(["bcftools", "view", "-h", "-Oz", "-o", str(empty), str(vcf)])
    run(["tabix", "-f", "-p", "vcf", str(empty)])
    info = {"vcf": str(vcf), "cadd": str(cadd), "empty": str(empty), "records": len(records)}
    (family / "done").write_text(json.dumps(info))
    return info


def checkout(work: Path, old_root: Path) -> Path:
    """A checkout with the old route's scripts (they find data/ by their real
    place) and a data/ folder holding the files the ACMG step names."""
    root = work / "checkout"
    if root.exists():
        return root
    shutil.copytree(old_root / "scripts", root / "scripts")
    for relative in ["alphamissense/alphamissense_mean_score.tsv", "repeats/repeats.bed.gz",
                     "MaveDB/assay_mapping_result.neutral_bands.tsv",
                     "ClinVar/clinvar_gene_stat.pkl",
                     "hpo/genes_to_phenotype.assertions.tsv.gz"]:
        path = root / "data" / relative
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text("stand-in\n")
    return root


def base_config(work: Path, family: dict) -> dict:
    resources = work / "resources"
    resources.mkdir(exist_ok=True)
    config = {}
    for key in RESOURCES:
        path = resources / key
        path.write_text("stand-in\n")
        config[key] = str(path)
    config.update({
        "assembly": "hg38", "threads_per_fam": 2, "af_cutoff": 0.01,
        "ped_file": str(FAMILY_PED), "clinvar_gene_stat": "clinvar_gene_stat.pkl",
        "repeat_region_file_name": "repeats.bed.gz", "alt_disease_vcf": family["vcf"],
        "extreme_rare_PAF": 0.0001, "exp_disease_incidence": None,
        "pp1_vcf": family["vcf"], "pp1_ped": str(FAMILY_PED),
        "apply_clingen_override": True, "pext_tissues": "Brain, Liver",
    })
    return config


def yaml_text(config: dict) -> str:
    lines = []
    for key, value in config.items():
        if value is None:
            lines.append(f"{key}: null")
        elif isinstance(value, bool):
            lines.append(f"{key}: {'true' if value else 'false'}")
        elif isinstance(value, (int, float)):
            lines.append(f"{key}: {value}")
        else:
            lines.append(f"{key}: {json.dumps(value)}")
    return "\n".join(lines) + "\n"


# Each case: the step, what result/ starts with, and the config changes.
def cases(family: dict) -> dict:
    table = "Chr\tStart\tRef\tAlt\tFeature\n22\t100\tA\tG\tENST1\n"
    done = "Chr\tStart\tRef\tAlt\tsplicing_lof\tsplicing_len_changing\tsplicing_canonical_offset\n"
    fam = "TEST_FAM"
    vcf_in = ("TEST.anno.vcf.gz", family["vcf"])
    cadd_in = ("TEST.anno.cadd.tsv", family["cadd"])
    fvcf = f"TEST.anno.{fam}.filtered.vcf.gz"
    ftab = f"TEST.anno.{fam}.filtered.tsv"
    return {
        "filtration-family": dict(step="filtration", family=fam, files=[vcf_in, cadd_in],
                                  vcf="TEST.anno.vcf.gz"),
        "filtration-no-ped": dict(step="filtration", files=[vcf_in, cadd_in],
                                  vcf="TEST.anno.vcf.gz", config={"ped_file": None}),
        "filtration-no-cadd": dict(step="filtration", family=fam, files=[vcf_in],
                                   vcf="TEST.anno.vcf.gz", config={"af_cutoff": 0.05}),
        "combine": dict(step="combine", files=[(fvcf, family["vcf"]),
                                               (f"TEST.anno.{fam}.filtered.cadd.tsv", family["cadd"])],
                        vcf=fvcf, cadd=f"TEST.anno.{fam}.filtered.cadd.tsv"),
        "combine-empty-vcf": dict(step="combine", files=[(fvcf, family["empty"]), cadd_in],
                                  vcf=fvcf, cadd="TEST.anno.cadd.tsv"),
        "combine-one-line-cadd": dict(step="combine", files=[(fvcf, family["vcf"]),
                                                             ("one.cadd.tsv", "## CADD\n")],
                                      vcf=fvcf, cadd="one.cadd.tsv"),
        "combine-default-hpo": dict(step="combine", files=[(fvcf, family["vcf"]), cadd_in],
                                    vcf=fvcf, cadd="TEST.anno.cadd.tsv",
                                    config={"hpo_assertions": None}),
        "splicing": dict(step="splicing", files=[(ftab, table)], table=ftab),
        "splicing-already-done": dict(step="splicing", files=[(ftab, done)], table=ftab),
        "splicing-empty": dict(step="splicing", files=[(ftab, "")], table=ftab),
        "utr": dict(step="utr", files=[(ftab, table)], table=ftab),
        "acmg-family-hg38": dict(step="acmg", family=fam, files=[(ftab, table)], table=ftab),
        "acmg-hg19-no-family": dict(step="acmg", files=[(ftab, table)], table=ftab, config={
            "assembly": "hg19", "apply_clingen_override": False, "pext_tissues": None,
            "exp_disease_incidence": 0.001, "cds_fasta_file": "/missing/cds.fa",
            "pp1_vcf": None, "relevant_gene_list": ""}),
        "acmg-gene-stat-path-missing": dict(step="acmg", family=fam, files=[(ftab, table)],
                                            table=ftab, config={"clinvar_gene_stat": "/abs/missing.pkl"}),
        "acmg-empty-values": dict(step="acmg", family=fam, files=[(ftab, table)], table=ftab,
                                  config={"pm1_regions_pkl": "", "alt_disease_vcf": "",
                                          "pp1_ped": "/missing.ped", "pext_tissues": ""}),
        "acmg-missing-resource": dict(step="acmg", family=fam, files=[(ftab, table)], table=ftab,
                                      config={"clinvar_aa_stat": "/missing/aa.pkl"}),
        "acmg-empty-table": dict(step="acmg", family=fam, files=[(ftab, "")], table=ftab),
        "prioritization": dict(step="prioritization", family=fam,
                               files=[(fvcf, family["vcf"]), cadd_in],
                               vcf=fvcf, cadd="TEST.anno.cadd.tsv"),
        "prioritization-empty-vcf": dict(step="prioritization", family=fam,
                                         files=[(fvcf, family["empty"]), cadd_in],
                                         vcf=fvcf, cadd="TEST.anno.cadd.tsv"),
    }


def stage(route: Path, case: dict, config: dict) -> Path:
    shutil.rmtree(route, ignore_errors=True)
    result = route / "result"
    result.mkdir(parents=True)
    for name, source in case["files"]:
        target = result / name
        if os.path.isfile(source):
            shutil.copyfile(source, target)
            if os.path.isfile(source + ".tbi"):
                shutil.copyfile(source + ".tbi", str(target) + ".tbi")
        else:
            target.write_text(source)
    config = dict(config, **case.get("config", {}))
    config["output_dir"] = str(result)
    config["input_vcf"] = str(route / "TEST.vcf.gz")
    path = route / "task.yaml"
    path.write_text(yaml_text(config))
    return path


def old_command(root: Path, case: dict, config: Path, result: Path) -> list:
    step, family = case["step"], case.get("family")
    if step == "filtration":
        command = ["bash", str(root / "scripts/filtration_vcf_per_fam.sh"), "-v", str(result / case["vcf"])]
        if family:
            command += ["-f", family]
        return command + ["-c", str(config)]
    if step == "prioritization":
        command = ["bash", str(root / "scripts/prioritization_vcf_per_fam.sh"), "-i", str(result / case["vcf"])]
        if family:
            command += ["-f", family]
        return command + ["-c", str(config), "-t", str(result / case["cadd"])]
    command = ["bash", str(FIXTURE), str(root), str(config), step]
    if step == "combine":
        return command + [str(result / case["vcf"]), str(result / case["cadd"])]
    command.append(str(result / case["table"]))
    if step == "acmg" and family:
        command.append(family)
    return command


def new_command(priva: str, root: Path, case: dict, config: Path, result: Path) -> list:
    command = [priva, "task", "step", "--step", case["step"], "--root", str(root),
               "--config", str(config)]
    for option in ("vcf", "cadd", "table"):
        if option in case:
            command += [f"--{option}", str(result / case[option])]
    if case.get("family"):
        command += ["--family", case["family"]]
    return command


def vcf_text(path: Path) -> str:
    text = run(["bcftools", "view", str(path)], capture_output=True, text=True).stdout
    return "\n".join(line for line in text.splitlines() if not line.startswith("##bcftools_"))


def outcome(route: Path, status: int) -> dict:
    result = route / "result"
    files = {}
    for path in sorted(result.rglob("*")):
        if not path.is_file():
            continue
        name = str(path.relative_to(result))
        if name.endswith((".tbi", ".csi")):
            files[name] = "index"
        elif name.endswith(".vcf.gz"):
            files[name] = vcf_text(path).replace(str(route), "ROUTE")
        else:
            files[name] = path.read_bytes()
    log = route / "priva.log"
    commands = log.read_text().replace(str(route), "ROUTE").splitlines() if log.exists() else []
    return {"status": status == 0, "commands": commands, "files": files}


def run_route(command: list, route: Path, priva_stand_in: Path, real_priva: str) -> int:
    temp = route / "tmp"
    temp.mkdir(exist_ok=True)
    env = dict(os.environ, PRIVA_BINARY=str(priva_stand_in), REAL_PRIVA=real_priva,
               FAKE_PRIVA_LOG=str(route / "priva.log"), TMPDIR=str(temp))
    with open(route / "stderr.txt", "w") as err:
        status = subprocess.run(command, cwd=route, env=env, stdout=err, stderr=err).returncode
    shutil.rmtree(temp, ignore_errors=True)
    return status


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    parser.add_argument("--priva", required=True)
    parser.add_argument("--work", required=True)
    parser.add_argument("--old-root", required=True)
    parser.add_argument("--case", action="append")
    args = parser.parse_args()
    old_root = Path(args.old_root).resolve()
    if not (old_root / "scripts" / "filtration_vcf_per_fam.sh").is_file():
        parser.error(f"--old-root {old_root} has no scripts/filtration_vcf_per_fam.sh (check out 841e200)")
    work = Path(args.work).resolve()
    work.mkdir(parents=True, exist_ok=True)
    priva = str(Path(args.priva).resolve())
    stand_in = work / "priva"
    stand_in.write_text(FAKE_PRIVA)
    stand_in.chmod(stand_in.stat().st_mode | stat.S_IXUSR | stat.S_IXGRP | stat.S_IXOTH)
    family = build_family(work)
    root = checkout(work, old_root)
    config = base_config(work, family)
    all_cases = cases(family)
    failed = 0
    for name in args.case or all_cases:
        case = all_cases[name]
        results = {}
        for route_name in ("old", "new"):
            route = work / name / route_name
            config_path = stage(route, case, config)
            result = route / "result"
            if route_name == "old":
                command = old_command(root, case, config_path, result)
            else:
                command = new_command(priva, root, case, config_path, result)
            status = run_route(command, route, stand_in, priva)
            results[route_name] = outcome(route, status)
        old, new = results["old"], results["new"]
        problems = []
        for key in ("status", "commands"):
            if old[key] != new[key]:
                problems.append(f"{key}:\n    old {old[key]}\n    new {new[key]}")
        if sorted(old["files"]) != sorted(new["files"]):
            problems.append(f"files: old {sorted(old['files'])} new {sorted(new['files'])}")
        for file in sorted(set(old["files"]) & set(new["files"])):
            if old["files"][file] != new["files"][file]:
                problems.append(f"{file} differs")
        kept = ""
        vcfs = [f for f in new["files"] if f.endswith(".vcf.gz") and new["files"][f] != "index"]
        for file in vcfs:
            records = sum(1 for line in new["files"][file].splitlines() if not line.startswith("#"))
            kept += f" {file}: {records} records;"
        state = "ok" if new["status"] else "refused"
        if problems:
            failed += 1
            print(f"DIFF {name} ({state})")
            for problem in problems:
                print("  " + problem)
        else:
            print(f"SAME {name} ({state}; {len(new['commands'])} priva commands;{kept}"
                  f" {len(new['files'])} files)")
    print(f"{failed} of {len(args.case or all_cases)} cases differ")
    return 1 if failed else 0


if __name__ == "__main__":
    sys.exit(main())
