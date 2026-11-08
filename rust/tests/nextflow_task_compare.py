#!/usr/bin/env python3
"""Compare `priva task` / `priva task-boundary` with the Python adapters.

Synthetic fixtures only (no genomic data). Each case runs both programs in
fresh folders on the same inputs and compares exit codes, refusal messages
and every file they write. Timestamps and paths are equal because both read
the same fixture files; per-run folders are mapped to one placeholder.

    python3 rust/tests/nextflow_task_compare.py --binary PRIVA --source CHECKOUT --work DIR
"""
import argparse
import gzip
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys

KEYS = ("clinvar_aa_stat", "clinvar_splice_stat", "clinvar_pm5_index", "all_intolerant_domains",
        "alphamissense_intolerant_domains", "pm1_regions_pkl", "alphamissense_vcf",
        "alphamissense_tranx_domain_map", "clinvar_patho_af_stat", "clinvar_patho_exon_af_stat",
        "interpro_mapping_pickle", "clingen_map", "gene_dosage_sensitivity",
        "hpo_condition_mechanism_json", "gene_nonlof_mechanism_json")
FAKE_FILTRATION = r"""#!/bin/bash
# Fake shell step: record what the task gave it, then write the outputs.
set -e
while [ $# -gt 0 ]; do case $1 in -v|-i) v=$2; shift 2;; -f) f=$2; shift 2;; -c) c=$2; shift 2;; -t) t=$2; shift 2;; *) shift;; esac; done
d=$(dirname "$v"); b=$(basename "$v" .vcf.gz)
yq -o=json . "$c" > "$d/../seen.$(basename "$0").json"
env | grep -E '^(PRIVA_BINARY|TMPDIR|TMP|TEMP|OMP_NUM_THREADS|OPENBLAS_NUM_THREADS|MKL_NUM_THREADS|NUMEXPR_NUM_THREADS|PYTHONDONTWRITEBYTECODE)=' | sort > "$d/../env.$(basename "$0").txt"
if [ "${0##*/}" = filtration_vcf_per_fam.sh ]; then
  cp "$v" "$d/$b.$f.filtered.vcf.gz"; cp "$v.tbi" "$d/$b.$f.filtered.vcf.gz.tbi"
  cp "$d/$b.cadd.tsv" "$d/$b.$f.filtered.cadd.tsv"
else
  printf 'a\n1\n' > "$d/$b.tsv"; printf 'b\n2\n' > "$d/$b.acmg.tsv"
fi
"""


def run(command, cwd, env=None):
    result = subprocess.run(command, cwd=cwd, capture_output=True, text=True, env=env)
    return result.returncode, result.stderr


def tree(root):
    """Every file below root with its bytes (directories and links by kind)."""
    found = {}
    for path in sorted(Path(root).rglob("*")):
        relative = str(path.relative_to(root))
        if path.name == "task.yaml":
            # Formatting differs by design; the fake steps record what yq reads.
            found[relative] = ("file", b"<compared through yq>")
        elif path.is_symlink():
            found[relative] = ("link", os.readlink(path))
        elif path.is_file():
            found[relative] = ("file", path.read_bytes())
        else:
            found[relative] = ("dir", None)
    return found


class Harness:
    def __init__(self, binary, source, work):
        self.binary, self.source, self.work = binary, Path(source), Path(work)
        self.python = [sys.executable, str(self.source / "bin/priva_task.py")]
        self.boundary = [sys.executable, str(self.source / "bin/priva_task_boundary.py")]
        self.failures, self.count = [], 0

    def record(self, name, ok, detail=""):
        self.count += 1
        print(("PASS " if ok else "FAIL ") + name + ("" if ok else ": " + detail[:2000]))
        if not ok:
            self.failures.append(name)

    def normalize(self, text, folders):
        for folder, label in folders:
            text = text.replace(str(folder), label)
        text = text.replace("priva_task.py manifest", "priva task manifest")
        # Usage errors: compare the argparse message, not the program name.
        return "\n".join(line.split(" error: ", 1)[-1] for line in text.splitlines())

    def both(self, name, args, setup=None, boundary=False, env=None, compare_files=True):
        """Run Python then Rust in their own folder; compare everything."""
        outcomes = []
        for label, program in (("python", self.boundary if boundary else self.python),
                               ("rust", [self.binary, "task-boundary" if boundary else "task"])):
            folder = self.work / name / label
            shutil.rmtree(folder, ignore_errors=True)
            folder.mkdir(parents=True)
            cwd = setup(folder) if setup else folder
            code, stderr = run(program + [a.replace("{folder}", str(folder)) for a in args], cwd, env)
            stderr = self.normalize(stderr, [(folder, "{folder}")])
            files = tree(folder) if compare_files else {}
            files = {k: (kind, self.normalize(v.decode("utf-8", "replace") if isinstance(v, bytes) else v,
                                              [(folder, "{folder}")]) if v is not None else v)
                     for k, (kind, v) in files.items()}
            outcomes.append((code, stderr.strip().splitlines()[-1:] if code else [], files))
        (pc, pe, pf), (rc, re_, rf) = outcomes
        detail = []
        if pc != rc or pe != re_:
            detail.append(f"exit/message python={pc} {pe} rust={rc} {re_}")
        for key in sorted(set(pf) | set(rf)):
            if pf.get(key) != rf.get(key):
                detail.append(f"file {key}: python={str(pf.get(key))[:600]!r} rust={str(rf.get(key))[:600]!r}")
        self.record(name, not detail, "\n".join(detail))
        return outcomes


def make_fixture(root, harness):
    """A fake legacy root, resources, a tiny synthetic VCF, PED and CADD."""
    root.mkdir(parents=True)
    source = root / "source"
    for name in ("main.nf", "nextflow.config"):
        (source / name).parent.mkdir(parents=True, exist_ok=True)
        (source / name).write_text("frozen\n")
    (source / "conf/cadd").mkdir(parents=True)
    (source / "conf/cadd.config").write_text("x\n")
    (source / "conf/cadd/profile.json").write_text("{}\n")
    (source / "bin").mkdir()
    (source / "bin/__pycache__").mkdir()
    (source / "bin/__pycache__/skip.pyc").write_bytes(b"\0")
    (source / "scripts").mkdir()
    for name in ("filtration_vcf_per_fam.sh", "prioritization_vcf_per_fam.sh"):
        (source / "scripts" / name).write_text(FAKE_FILTRATION)
    for relative in ("data/hpo/genes_to_phenotype.assertions.tsv.gz",
                     "data/alphamissense/alphamissense_mean_score.tsv", "data/ClinVar/gene.pkl",
                     "data/repeats/repeats.bed.gz", "data/hgnc/non_alt_loci_set.tsv",
                     "data/dispensable_genes/dispensable_gene_list.txt",
                     "data/MaveDB/assay_mapping_result.neutral_bands.tsv"):
        (source / relative).parent.mkdir(parents=True, exist_ok=True)
        (source / relative).write_bytes(b"synthetic")
    resources = root / "resources"
    (resources / "tree/sub").mkdir(parents=True)
    (resources / "tree/sub/a.tsv").write_text("a\n")
    (resources / "tree/b.tsv").write_text("b\n")
    (resources / "tree/.hidden").write_text("h\n")
    os.symlink(resources / "tree/b.tsv", resources / "tree/link.tsv")
    (resources / "one.pkl").write_bytes(b"resource")
    (resources / "scores.vcf.gz").write_bytes(b"scores")
    (resources / "scores.vcf.gz.tbi").write_bytes(b"index")
    (resources / "controls.vcf.gz").write_bytes(b"controls")
    (resources / "controls.vcf.gz.tbi").write_bytes(b"cindex")
    (resources / "coding.fa.gz").write_bytes(b"fasta")
    (resources / "coding.cds_cache.pkl.gz").write_bytes(b"cache")
    priva = root / "priva"
    os.symlink(harness.binary, priva)
    config = {key: str(resources / "one.pkl") for key in KEYS}
    config.update(alphamissense_vcf=str(resources / "scores.vcf.gz"),
                  vep_like_tree=str(resources / "tree"), interpro_mapping_pickle=str(resources / "tree"),
                  assembly="GRCh38", clinvar_gene_stat="gene.pkl", repeat_region_file_name="repeats.bed.gz",
                  priva_binary=str(priva), threads=2, threads_per_fam=1, af_cutoff=0.05,
                  apply_clingen_override=False, control_vcf=str(resources / "controls.vcf.gz"),
                  cds_fasta_file=str(resources / "coding.fa.gz"), hub_vcf_file="/shared/hub.vcf.gz",
                  note="text with 'quotes' and é")
    import yaml
    (root / "legacy.yaml").write_text(yaml.safe_dump(config, sort_keys=False)
                                      + "flag_yes: yes  # PyYAML true\n")
    (root / "family.ped").write_text("OTHER\tother\t0\t0\t2\t2\nFAM\tchild\tdad\tmom\t1\t2\n"
                                     "FAM\tmom\t0\t0\t2\t1\nFAM\tdad\t0\t0\t1\t1\n")
    body = ("##fileformat=VCFv4.2\n##contig=<ID=chr1,length=1000>\n"
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tchild\tdad\tmom\tother\n"
            "chr1\t10\t.\tA\tG\t50\tPASS\t.\tGT\t0/1\t0/0\t0/1\t0/0\n")
    (root / "input.vcf").write_text(body)
    subprocess.run(["bgzip", "-f", str(root / "input.vcf")], check=True)
    subprocess.run(["tabix", "-p", "vcf", str(root / "input.vcf.gz")], check=True)
    (root / "input.cadd.tsv").write_text("#Chrom\tPos\tRef\tAlt\tRawScore\tPHRED\tFeatureID\n1\t10\tA\tG\t1\t10\tX\n")
    return source


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--binary", required=True)
    parser.add_argument("--source", required=True)
    parser.add_argument("--work", required=True)
    args = parser.parse_args()
    harness = Harness(str(Path(args.binary).resolve()), args.source, args.work)
    work = Path(args.work)
    shutil.rmtree(work, ignore_errors=True)
    fixture = work / "fixture"
    source = make_fixture(fixture, harness)
    f = lambda name: str(fixture / name)
    common = ["--config", f("legacy.yaml"), "--legacy-root", str(source)]

    # 1. Manifests: identical bytes, in both identity modes.
    for mode, extra in (("annotated", []), ("annotated", ["--checksum"]),
                        ("annotated", ["--cadd-hub-mode", "off"])):
        name = "manifest_" + mode + "".join(extra).replace("-", "_")
        harness.both(name, ["manifest"] + common + ["--mode", mode, "--output", "{folder}/m.json"] + extra)
    manifest = work / "shared_manifest.json"
    subprocess.run(harness.python + ["manifest"] + common + ["--mode", "annotated", "--output", str(manifest)], check=True)
    digest = subprocess.run([sys.executable, "-c", "import sys,importlib.util as u;"
                             "s=u.spec_from_file_location('t',sys.argv[1]);m=u.module_from_spec(s);"
                             "s.loader.exec_module(m);print(m.source_digest(sys.argv[2]))",
                             str(harness.source / "bin/priva_task.py"), str(source)],
                            capture_output=True, text=True, check=True).stdout.strip()

    # 2. Launch markers and refusals.
    launch = ["launch", "--vcf", f("input.vcf.gz"), "--vcf-index", f("input.vcf.gz.tbi"),
              "--ped", f("family.ped"), "--config", f("legacy.yaml"), "--resources", str(manifest),
              "--family", "FAM", "--assembly", "hg38", "--mode", "annotated", "--legacy-root", str(source),
              "--source-digest", digest, "--software-id", "fixture env", "--vcf-hub-mode", "off",
              "--cadd-hub-mode", "copy", "--annotation-cpus", "1", "--family-cpus", "2",
              "--cadd", f("input.cadd.tsv")]
    harness.both("launch_fresh", launch + ["--outdir", "{folder}/out"])
    harness.both("launch_stub", launch + ["--outdir", "{folder}/out", "--stub"])
    def nonempty(folder):
        (folder / "out").mkdir()
        (folder / "out/other").write_text("x")
        return folder
    harness.both("launch_nonempty_outdir", launch + ["--outdir", "{folder}/out"], setup=nonempty)
    def foreign(folder):
        (folder / "out").mkdir()
        (folder / "out/.priva-run.json").write_text('{"run_identity": "other"}\n')
        return folder
    harness.both("launch_foreign_identity", launch + ["--outdir", "{folder}/out"], setup=foreign)
    def same(folder):
        subprocess.run(harness.python + launch + ["--outdir", str(folder / "out")], check=True)
        return folder
    harness.both("launch_identical_resume", launch + ["--outdir", "{folder}/out"], setup=same)
    swap = lambda items, key, value: [value if i and items[i - 1] == key else x for i, x in enumerate(items)]
    for name, changed in (("launch_bad_digest", swap(launch, "--source-digest", "0" * 64)),
                          ("launch_wrong_assembly", swap(launch, "--assembly", "hg19")),
                          ("launch_raw_gated", swap(launch, "--mode", "raw")),
                          ("launch_too_few_family_cpus", swap(launch, "--family-cpus", "0")),
                          ("launch_bad_int", swap(launch, "--family-cpus", "two"))):
        harness.both(name, changed + ["--outdir", "{folder}/out"])

    # 3. Validation bundles (real tools when present, and stub).
    validate = ["validate", "--vcf", f("input.vcf.gz"), "--vcf-index", f("input.vcf.gz.tbi"),
                "--ped", f("family.ped"), "--config", f("legacy.yaml"), "--resources", str(manifest),
                "--family", "FAM", "--base", "input", "--assembly", "hg38", "--mode", "annotated",
                "--legacy-root", str(source), "--source-digest", digest, "--software-id", "fixture",
                "--vcf-hub-mode", "off", "--cadd-hub-mode", "copy", "--cadd", f("input.cadd.tsv")]
    harness.both("validate_stub", validate + ["--stub"])
    outcomes = harness.both("validate_real", validate)
    for name, changed in (("validate_unsafe_family", swap(validate, "--family", "F;AM")),
                          ("validate_absent_family", swap(validate, "--family", "NONE")),
                          ("validate_boundary_inactive", validate + ["--task-boundary", "bwrap"]),
                          ("validate_bad_cadd", swap(validate, "--cadd", f("family.ped"))),
                          ("validate_missing_index", swap(validate, "--vcf-index", f("absent.tbi")))):
        harness.both(name, changed)

    # 4. Task phases on a Python-made bundle: stub, preannotated import and
    #    the two family phases with fake shell steps reading the task config.
    bundle = work / "bundle_real"
    bundle.mkdir()
    subprocess.run(harness.python + validate, cwd=bundle, check=True)
    stub_bundle = work / "bundle_stub"
    stub_bundle.mkdir()
    subprocess.run(harness.python + validate + ["--stub"], cwd=stub_bundle, check=True)
    phase_inputs = {"annotation": [], "filtration": ["--cadd", f("input.cadd.tsv")],
                    "prioritization": ["--cadd", f("input.cadd.tsv")]}
    for phase, extra in phase_inputs.items():
        for label, chosen, flags in (("stub", stub_bundle, ["--stub"]), ("real", bundle, [])):
            if label == "real" and phase == "annotation":
                extra = ["--cadd", f("input.cadd.tsv")]
            harness.both(f"run_{phase}_{label}", ["run", "--phase", phase, "--bundle", str(chosen / "bundle"),
                         "--vcf", f("input.vcf.gz"), "--vcf-index", f("input.vcf.gz.tbi"), "--cpus", "2"]
                         + extra + flags)
    harness.both("run_mixed_stub", ["run", "--phase", "filtration", "--bundle", str(bundle / "bundle"),
                 "--vcf", f("input.vcf.gz"), "--vcf-index", f("input.vcf.gz.tbi"), "--cpus", "1",
                 "--cadd", f("input.cadd.tsv"), "--stub"])
    harness.both("run_missing_cadd", ["run", "--phase", "filtration", "--bundle", str(bundle / "bundle"),
                 "--vcf", f("input.vcf.gz"), "--vcf-index", f("input.vcf.gz.tbi"), "--cpus", "1"])
    harness.both("run_too_few_cpus", ["run", "--phase", "prioritization", "--bundle", str(bundle / "bundle"),
                 "--vcf", f("input.vcf.gz"), "--vcf-index", f("input.vcf.gz.tbi"), "--cpus", "0",
                 "--cadd", f("input.cadd.tsv")])
    # Resource changed after the bundle was made.
    (fixture / "resources/tree/sub/a.tsv").write_text("changed size\n")
    harness.both("run_resource_changed", ["run", "--phase", "filtration", "--bundle", str(bundle / "bundle"),
                 "--vcf", f("input.vcf.gz"), "--vcf-index", f("input.vcf.gz.tbi"), "--cpus", "1",
                 "--cadd", f("input.cadd.tsv")])

    # 5. The task boundary check.
    def task(folder, tweak=None):
        (folder / "task").mkdir()
        for name in (".command.sh", ".command.run"):
            (folder / "task" / name).write_text("true\n")
        if tweak:
            tweak(folder / "task", folder)
        return folder / "task"
    for name, tweak in (("ok", None),
                        ("hardlink", lambda t, f: (f / "c").write_text("c") or os.link(f / "c", t / "alias")),
                        ("fifo", lambda t, f: os.mkfifo(t / "fifo")),
                        ("missing_control", lambda t, f: (t / ".command.run").unlink()),
                        ("staged_dir_link", lambda t, f: (f / "o").mkdir() or os.mkfifo(f / "o/p")
                         or os.symlink(f / "o", t / "in")),
                        ("scratch_link", lambda t, f: os.symlink(f, t / ".priva_task_tmp"))):
        harness.both(f"boundary_check_{name}", ["check", "--task-dir", "."],
                     setup=lambda folder, tw=tweak: task(folder, tw), boundary=True)
    harness.both("boundary_wrong_cwd", ["check", "--task-dir", "{folder}/task"],
                 setup=lambda folder: (task(folder), folder)[1], boundary=True)
    harness.both("boundary_exec_empty", ["exec", "--task-dir", "."],
                 setup=lambda folder: task(folder), boundary=True)
    print(f"{harness.count - len(harness.failures)} of {harness.count} cases match")
    return 1 if harness.failures else 0


if __name__ == "__main__":
    sys.exit(main())
