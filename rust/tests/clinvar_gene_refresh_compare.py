#!/usr/bin/env python3
"""Compare `priva clinvar-gene-refresh` with its unchanged Python source.

check --binary PRIVA --repo CHECKOUT --work-root NEW_DIRECTORY [--network] [--only CASE ...]
oracle --producer SCRIPT --exception JSON -- ARGUMENTS (internal reference launcher)

Each case builds one fixture folder: a source root holding copies of the
scientific scripts the profile binds, a recording install_utils.sh with the
four ClinVar stages, the resource registry and a small bundled baseline; a
managed cache with generated VCFs (pysam BGZF and tabix); and a configuration.
Every fixture file gets a fixed mtime. The Python arm runs
CHECKOUT/scripts/clinvar_gene_refresh.py through the oracle launcher, which
records the class and text of the exception main() reports; the folder is then
moved aside and built again at the same path for the native arm, so absolute
paths, sizes, digests and mtimes agree. Both arms run `python` from PATH, the
same interpreter the native profile observes.

Compared exactly: exit status, stdout, the exception category and message, and
the folder after the run (every entry's kind, mode and content). UTC
timestamps, retained-generation folder names and the device, inode and ctime
fields of retirement identities are written as placeholders. Argument errors
exit 2 in both. Cases marked network reach NCBI (a .md5 manifest and one 1 MiB
range of clinvar.vcf.gz each); they run only with --network. Neither arm
acquires full NCBI files, so acquire and a rebuilding install are untested.
Writes WORK_ROOT/complete.json and exits 1 on any mismatch.
"""
from __future__ import annotations

import argparse
import hashlib
import importlib.util
import json
import os
from pathlib import Path
import re
import shutil
import subprocess
import sys
import time

sys.dont_write_bytecode = True
PRODUCER_SHA256 = "ee238e18dfd582b50b62fad346aa89557c88200a52f06198ae00be12c8156d86"
MTIME = 1_700_000_000_123_456_789
OLD = 1_600_000_000_000_000_000
# A recent bundled baseline: its three-month boundary lies after the latest ClinVar release.
RECENT = time.time_ns()
SCIENTIFIC = ("clinvar_stat_variants.py", "combine_annotations.py", "annotation_vcf.sh",
              "clinvar_annotation_contract.py", "common_bash_utils.sh", "stat_aachange_clinvar.py",
              "hpo_penetrance.py", "clinvar_gene_refresh.py")
TIME = re.compile(r"\d{4}-\d{2}-\d{2}T\d{2}:\d{2}:\d{2}(?:\.\d{6})?\+00:00")
ARCHIVE = re.compile(r"\d{8}T\d{6}Z-[A-Za-z0-9_]+")
RAW = re.compile(r'"(dev|ino|ctime_ns)": "\d+"')
VEP116 = '##VEP="v116" API="v116" assembly="GRCh38.p14" cache="x"'
INSTALLER = r'''#!/bin/bash
# Recording stand-in for the ClinVar stages of install_utils.sh.
record() { printf '%s\n' "$*" >> "${FIXTURE_LOG:?}"; }
function ClinVar_VCF_deploy() {
    record "ClinVar_VCF_deploy $*"
}


function _ClinVar_VCF_deploy() {
    record "_ClinVar_VCF_deploy $*"
    [[ -z ${FIXTURE_DEPLOY_STATUS:-} ]] || return "$FIXTURE_DEPLOY_STATUS"
    if [[ -n ${FIXTURE_DEPLOY_OUTPUT:-} ]]; then
        cp "$FIXTURE_DEPLOY_SOURCE" "$FIXTURE_DEPLOY_OUTPUT"
        touch -d @1700000000 "$FIXTURE_DEPLOY_OUTPUT"
    fi
}


function ClinVar_AA_stat () {
    record "ClinVar_AA_stat $*"
}


function ClinVar_Gene_stat_regenerate () {
    record "ClinVar_Gene_stat_regenerate $*"
}


function ClinVar_Gene_stat () {
    record "ClinVar_Gene_stat $*"
}

"$@"
'''


def load(path, name):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def write_json(path, value):
    Path(path).write_text(json.dumps(value, sort_keys=True, indent=2) + "\n")


def oracle(args):
    """Run the unchanged main() and record the exception it reports."""
    sys.path.insert(0, str(Path(args.producer).parent))
    module = load(args.producer, "original_clinvar_gene_refresh")
    record = {}
    original_exit = module.argparse.ArgumentParser.exit

    def observed_exit(parser, status=0, message=None):
        kind, value, _ = sys.exc_info()
        if status == 1 and kind is not None:
            record.update(category=kind.__name__, message=str(value))
            write_json(args.exception, record)
        return original_exit(parser, status, message)

    module.argparse.ArgumentParser.exit = observed_exit
    sys.argv = [str(args.producer), *args.arguments]
    try:
        module.main()
    except SystemExit:
        raise
    except BaseException as exc:
        record.update(category=type(exc).__name__, message=str(exc), uncaught=True)
        write_json(args.exception, record)
        raise
    raise SystemExit(0)


# ---------------------------------------------------------------- fixtures

class Fixture:
    """One case folder; every helper returns the path it wrote."""

    def __init__(self, root, repo):
        self.root, self.repo = root, repo
        self.source = root / "source"
        self.cache = root / "cache"
        self.clinvar = self.cache / "clinvar"
        self.target = self.clinvar / "clinvar.hg38.vep.vcf.gz"
        self.config = root / "config.yaml"
        self.log = root / "stages.log"
        self.env = {"FIXTURE_LOG": str(self.log)}
        self.times = {}
        scripts = self.source / "scripts"
        scripts.mkdir(parents=True)
        for name in SCIENTIFIC:
            shutil.copyfile(repo / "scripts" / name, scripts / name)
        (scripts / "install_utils.sh").write_text(INSTALLER)
        shutil.copyfile(repo / "cache_resources.json", self.source / "cache_resources.json")
        self.baseline = self.write(self.source / "data/ClinVar/clinvar_2star_stats.pkl.gz", b"baseline")
        self.clinvar.mkdir(parents=True)
        self.settings = {"assembly": "hg38", "vep_installed_version": 116, "cache_root": str(self.cache),
                         "clinvar_vcf_dir": str(self.clinvar), "threads": 1}

    def write(self, path, data):
        path = Path(path)
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_bytes(data if isinstance(data, bytes) else data.encode())
        return path

    def vcf(self, path, *, header=VEP116, annotated=False, date="20261027", index=True, plain=False):
        import pysam
        sys.path.insert(0, str(self.repo / "scripts"))
        import clinvar_annotation_contract as contract
        text = f"##fileformat=VCFv4.2\n##fileDate={date}\n##contig=<ID=chr1,length=1000>\n"
        if header:
            text += header + "\n"
        if annotated:
            fields = sorted(contract.AA_CSQ_FIELDS | contract.GENE_CSQ_FIELDS)
            text += '##INFO=<ID=CSQ,Number=.,Type=String,Description="Format: ' + "|".join(fields) + '">\n'
            for name in ("CLNSIG", "CLNREVSTAT", "ALLELEID"):
                text += f'##INFO=<ID={name},Number=1,Type=String,Description="fixture">\n'
        text += "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\nchr1\t10\t.\tA\tC\t.\t.\t.\n"
        path = Path(path)
        path.parent.mkdir(parents=True, exist_ok=True)
        if plain:
            path.write_text(text)
            return path
        with pysam.BGZFile(str(path), "wb") as output:
            output.write(text.encode())
        if index:
            pysam.tabix_index(str(path), preset="vcf", force=True)
        return path

    def registry(self, change):
        path = self.source / "cache_resources.json"
        value = json.loads(path.read_text())
        change(value)
        path.write_text(json.dumps(value, indent=2))

    def save(self, **settings):
        import yaml
        self.settings.update(settings)
        self.config.write_text(yaml.safe_dump({k: v for k, v in self.settings.items() if v is not ...},
                                              sort_keys=False))

    def profile_inputs(self):
        """Selected resources the scientific profile records."""
        cache = self.cache
        self.write(cache / "reference.fa", ">chr1\nACGT\n")
        self.write(cache / "reference.fa.fai", "chr1\t4\t6\t4\t5\n")
        for name in ("UTRAnnotator", "SpliceAI", "NMD"):
            self.write(cache / f"plugins/{name}.pm", f"package {name};\n1;\n")
        for name in ("LoF.pm", "utr_splice.pl", "README.md", ".hidden.pl"):
            self.write(cache / f"loftee/{name}", name)
        self.write(cache / "gnomad/gnomad.chrX.vcf.gz", "x")
        self.write(cache / "gnomad/gnomad.chr1.vcf.gz", "1")
        self.write(cache / "gnomad/gnomad.chr1.vcf.gz.tbi", "1i")
        self.write(cache / "vep/homo_sapiens_merged/116_GRCh38/info.txt", "info\n")
        self.write(cache / "pext/tissues/Liver.bw", "liver")
        self.write(self.source / "data/loeuf/loeuf.tsv.gz", "loeuf")
        self.write(cache / "big.bin", b"\0" * (1024 * 1024 + 1))
        self.save(ref_genome=str(cache / "reference.fa"), vep_plugins_dir=str(cache / "plugins"),
                  vep_cache_dir=str(cache / "vep"), loftee_repo=str(cache / "loftee"),
                  gnomad_vcf_chrX=str(cache / "gnomad/gnomad.chrX.vcf.gz"), loeuf_prescore="loeuf.tsv.gz",
                  conservation_file=str(cache / "big.bin"), pext_enabled=True,
                  pext_tissue_bw_dir=str(cache / "pext"), pext_tissues="Liver, Brain",
                  human_ancestor_fasta=str(cache / "missing.fa"))

    def spliceai(self):
        snv = self.vcf(self.cache / "spliceai/snv.vcf.gz", header=None)
        indel = self.vcf(self.cache / "spliceai/indel.vcf.gz", header=None)
        self.save(spliceai_snv_prescore=str(snv), spliceai_indel_prescore=str(indel))

    def finish(self):
        for folder, _, names in os.walk(self.root):
            for name in names:
                path = Path(folder) / name
                if not path.is_symlink():
                    os.utime(path, ns=(MTIME, self.times.get(path, MTIME)))


def common(f):
    return ["--config", str(f.config), "--source-root", str(f.source)]


def catalogue():
    """Case name -> (build(fixture) returning the argument list, network)."""
    cases = {}

    def case(name, network=False):
        def register(build):
            cases[name] = (build, network)
            return build
        return register

    # ---- file-date
    @case("file_date_gzip")
    def _(f):
        return ["file-date", "--input", str(f.vcf(f.root / "a.vcf.gz"))]

    @case("file_date_plain_dashed")
    def _(f):
        return ["file-date", "--input", str(f.vcf(f.root / "a.vcf", date="2026-10-27", plain=True))]

    for label, date in (("conflict", "20261027\n##fileDate=20261028"), ("bad_form", "2026/10/27"),
                        ("impossible_day", "2026-02-30"), ("month_13", "2026-13-01"),
                        ("compact_month_13", "20261301"), ("day_32", "2026-10-32"), ("year_0", "0000-01-01")):
        @case("file_date_" + label)
        def _(f, date=date):
            return ["file-date", "--input", str(f.vcf(f.root / "a.vcf", date=date, plain=True))]

    @case("file_date_missing_date")
    def _(f):
        f.write(f.root / "a.vcf", "##fileformat=VCFv4.2\n#CHROM\tPOS\n")
        return ["file-date", "--input", str(f.root / "a.vcf")]

    @case("file_date_records_first")
    def _(f):
        f.write(f.root / "a.vcf", "chr1\t1\n")
        return ["file-date", "--input", str(f.root / "a.vcf")]

    @case("file_date_empty")
    def _(f):
        f.write(f.root / "a.vcf", "")
        return ["file-date", "--input", str(f.root / "a.vcf")]

    @case("file_date_absent")
    def _(f):
        return ["file-date", "--input", str(f.root / "absent.vcf")]

    @case("file_date_bad_gzip")
    def _(f):
        f.write(f.root / "a.vcf.gz", b"\x1f\x8b\x08\x00garbage")
        return ["file-date", "--input", str(f.root / "a.vcf.gz")]

    @case("file_date_truncated_gzip")
    def _(f):
        data = f.vcf(f.root / "whole.vcf.gz", index=False).read_bytes()
        f.write(f.root / "a.vcf.gz", data[:40])
        return ["file-date", "--input", str(f.root / "a.vcf.gz")]

    # ---- candidate-route
    @case("route")
    def _(f):
        return ["candidate-route", "--source-root", str(f.source)]

    @case("route_missing_entry")
    def _(f):
        path = f.source / "scripts/annotation_vcf.sh"
        path.write_text(path.read_text().replace("function anno_VEP_cache_candidate() {", "function renamed() {"))
        return ["candidate-route", "--source-root", str(f.source)]

    # ---- plan without NCBI
    @case("plan_explicit_preserved")
    def _(f):
        stats = f.write(f.cache / "stats/selected.rkyv", "stats")
        f.save(clinvar_gene_stat=str(stats))
        return ["plan", *common(f)]

    @case("plan_explicit_output")
    def _(f):
        stats = f.write(f.cache / "stats/selected.rkyv", "stats")
        f.save(clinvar_gene_stat=str(stats))
        return ["plan", *common(f), "--output", str(f.root / "plans/plan.json")]

    @case("plan_explicit_empty")
    def _(f):
        f.save(clinvar_gene_stat=str(f.write(f.cache / "stats/selected.rkyv", "")))
        return ["plan", *common(f)]

    @case("plan_explicit_absent")
    def _(f):
        f.save(clinvar_gene_stat=str(f.cache / "stats/absent.rkyv"))
        return ["plan", *common(f)]

    @case("plan_disabled")
    def _(f):
        f.save()
        f.registry(lambda v: v["resources"]["clinvar_gene_statistics"]["refresh"].update(automatic=False))
        return ["plan", *common(f)]

    @case("plan_bad_assembly")
    def _(f):
        f.save(assembly="hg18")
        return ["plan", *common(f)]

    @case("plan_bad_major")
    def _(f):
        f.save(vep_installed_version="v116")
        return ["plan", *common(f)]

    @case("plan_not_mapping")
    def _(f):
        f.config.write_text("- a\n- b\n")
        return ["plan", *common(f)]

    @case("plan_duplicate_key")
    def _(f):
        f.config.write_text("assembly: hg38\nassembly: hg19\n")
        return ["plan", *common(f)]

    registry_changes = {
        "missing": None,
        "not_json": "text",
        "schema": lambda v: v.update(schema_version=2),
        "no_refresh": lambda v: v["resources"]["clinvar_gene_statistics"].pop("refresh"),
        "automatic_text": lambda v: v["resources"]["clinvar_gene_statistics"]["refresh"].update(automatic="yes"),
        "months_zero": lambda v: v["resources"]["clinvar_gene_statistics"]["refresh"].update(threshold_months=0),
        "months_bool": lambda v: v["resources"]["clinvar_gene_statistics"]["refresh"].update(threshold_months=True),
        "template_two": lambda v: v["resources"]["clinvar_vcf"]["source"].update(
            latest_url_template="https://ftp.ncbi.nlm.nih.gov/{x}/vcf_{assembly}/clinvar.vcf.gz"),
        "template_name": lambda v: v["resources"]["clinvar_vcf"]["source"].update(
            latest_url_template="https://ftp.ncbi.nlm.nih.gov/vcf_{assembly}/clinvar.vcf"),
        "not_ncbi": lambda v: v["resources"]["mane_select"]["source"].update(url="http://example.org/x"),
        "resources_list": lambda v: v.update(resources=[]),
    }
    for label, change in registry_changes.items():
        @case("plan_registry_" + label)
        def _(f, change=change):
            f.save()
            path = f.source / "cache_resources.json"
            if change is None:
                path.unlink()
            elif isinstance(change, str):
                path.write_text(change)
            else:
                f.registry(change)
            return ["plan", *common(f)]

    # ---- plan with NCBI
    @case("plan_rebuild_profile", network=True)
    def _(f):
        f.profile_inputs()
        f.times[f.baseline] = OLD
        return ["plan", *common(f)]

    @case("plan_before_boundary", network=True)
    def _(f):
        f.save()
        f.times[f.baseline] = RECENT
        return ["plan", *common(f)]

    @case("plan_missing_baseline", network=True)
    def _(f):
        f.save(priva_python_hash_seed=7, pext_enabled=False)
        f.baseline.unlink()
        return ["plan", *common(f)]

    @case("plan_http_error", network=True)
    def _(f):
        f.save()
        f.times[f.baseline] = OLD
        f.registry(lambda v: v["resources"]["clinvar_vcf"]["source"].update(
            latest_url_template="https://ftp.ncbi.nlm.nih.gov/pub/clinvar/absent_{assembly}/clinvar.vcf.gz"))
        return ["plan", *common(f)]

    # ---- candidate begin/complete
    def candidate(f, annotated=False):
        f.profile_inputs()
        f.spliceai()
        vcf = f.vcf(f.target, annotated=annotated)
        f.save(clinvar_vcf=str(vcf))
        return vcf

    @case("candidate_begin")
    def _(f):
        vcf = candidate(f)
        return ["candidate-begin", *common(f), "--input-vcf", str(vcf)]

    @case("candidate_begin_outside_cache")
    def _(f):
        candidate(f)
        other = f.vcf(f.root / "other.vcf.gz")
        return ["candidate-begin", *common(f), "--input-vcf", str(other)]

    @case("candidate_begin_relative_root")
    def _(f):
        vcf = candidate(f)
        f.save(cache_root="cache")
        return ["candidate-begin", *common(f), "--input-vcf", str(vcf)]

    @case("candidate_begin_no_spliceai")
    def _(f):
        vcf = candidate(f)
        f.save(spliceai_indel_prescore=None)
        return ["candidate-begin", *common(f), "--input-vcf", str(vcf)]

    @case("candidate_complete")
    def _(f):
        vcf = candidate(f, annotated=True)
        return [("candidate-begin", *common(f), "--input-vcf", str(vcf)),
                ["candidate-complete", *common(f), "--input-vcf", str(vcf)]]

    @case("candidate_complete_changed")
    def _(f):
        vcf = candidate(f, annotated=True)
        return [("candidate-begin", *common(f), "--input-vcf", str(vcf)),
                ("append", str(f.config), "threads: 2\n"),
                ["candidate-complete", *common(f), "--input-vcf", str(vcf)]]

    @case("candidate_complete_without_begin")
    def _(f):
        vcf = candidate(f, annotated=True)
        return ["candidate-complete", *common(f), "--input-vcf", str(vcf)]

    @case("candidate_complete_unannotated")
    def _(f):
        vcf = candidate(f)
        return [("candidate-begin", *common(f), "--input-vcf", str(vcf)),
                ["candidate-complete", *common(f), "--input-vcf", str(vcf)]]

    # ---- check-header and deployment
    def owned(f, header=VEP116):
        f.save()
        return f.vcf(f.target, header=header)

    @case("check_header_owned")
    def _(f):
        owned(f)
        return ["check-header", *common(f)]

    for label, header in (("old_release", '##VEP="v115" assembly="GRCh38"'),
                          ("hg19", '##VEP="v116" assembly="GRCh37.p13"'),
                          ("missing", None),
                          ("two_vep", VEP116 + "\n" + VEP116),
                          ("unclosed_quote", '##VEP="v116 assembly=GRCh38')):
        @case("check_header_" + label)
        def _(f, header=header):
            owned(f, header)
            return ["check-header", *common(f)]

    @case("check_header_absent")
    def _(f):
        f.save()
        return ["check-header", *common(f)]

    @case("check_header_linked")
    def _(f):
        real = f.vcf(f.root / "elsewhere.vcf.gz")
        f.save()
        f.target.symlink_to(real)
        return ["check-header", *common(f)]

    @case("check_header_owned_outside_root")
    def _(f):
        owned(f)
        f.save(cache_root=str(f.root / "other"))
        return ["check-header", *common(f)]

    def deploy_env(f, source=None):
        f.env.update(FIXTURE_DEPLOY_OUTPUT=str(f.target),
                     FIXTURE_DEPLOY_SOURCE=str(source or f.vcf(f.root / "fresh.vcf.gz")))

    @case("deploy_owned_current")
    def _(f):
        owned(f)
        deploy_env(f)
        return ["deploy", *common(f)]

    @case("deploy_owned_absent")
    def _(f):
        f.save()
        deploy_env(f)
        return ["deploy", *common(f)]

    @case("deploy_owned_obsolete_retained")
    def _(f):
        owned(f, '##VEP="v115" assembly="GRCh38"')
        f.write(f.clinvar / "clinvar.GRCh38.vcf.gz", "raw")
        f.write(f.clinvar / "clinvar.hg38.vep.gnomAD.log", "log")
        f.write(f.clinvar / "clinvar.hg38.vep.vcf", "plain")
        f.times[f.clinvar / "clinvar.hg38.vep.vcf"] = OLD
        deploy_env(f)
        return ["deploy", *common(f)]

    @case("deploy_owned_newer_plain_sibling")
    def _(f):
        owned(f)
        f.vcf(f.clinvar / "clinvar.hg38.vep.vcf", header='##VEP="v110" assembly="GRCh38"', plain=True)
        f.times[f.clinvar / "clinvar.hg38.vep.vcf"] = MTIME + 1000
        deploy_env(f)
        return ["deploy", *common(f)]

    @case("deploy_owned_pinned_source")
    def _(f):
        owned(f, None)
        pinned = f.vcf(f.cache / "pinned/clinvar_raw.vcf.gz", header=None)
        deploy_env(f)
        return ["deploy", *common(f), "--source-vcf", str(pinned)]

    @case("deploy_pinned_source_alias")
    def _(f):
        owned(f, None)
        pinned = f.vcf(f.clinvar / "clinvar.hg38.raw.vcf.gz", header=None)
        deploy_env(f)
        return ["deploy", *common(f), "--source-vcf", str(pinned)]

    @case("deploy_stage_fails")
    def _(f):
        f.save()
        f.env.update(FIXTURE_DEPLOY_STATUS="3")
        return ["deploy", *common(f)]

    @case("deploy_output_obsolete")
    def _(f):
        f.save()
        deploy_env(f, f.vcf(f.root / "fresh.vcf.gz", header='##VEP="v115" assembly="GRCh38"'))
        return ["deploy", *common(f)]

    @case("deploy_explicit_current")
    def _(f):
        explicit = f.vcf(f.cache / "explicit/clinvar.custom.vcf.gz")
        f.save(clinvar_vcf=str(explicit))
        return ["deploy", *common(f)]

    @case("deploy_explicit_incompatible")
    def _(f):
        explicit = f.vcf(f.cache / "explicit/clinvar.custom.vcf.gz", header='##VEP="v99" assembly="GRCh38"')
        f.save(clinvar_vcf=str(explicit))
        return ["deploy", *common(f)]

    @case("deploy_owned_earlier_vcf_component")
    def _(f):
        folder = f.cache / "a.vcf.d"
        folder.mkdir(parents=True)
        f.save(clinvar_vcf_dir=str(folder))
        return ["deploy", *common(f)]

    # ---- retire-owned and retire-shards
    @case("retire_owned")
    def _(f):
        owned(f, '##VEP="v115" assembly="GRCh38"')
        return ["retire-owned", *common(f)]

    @case("retire_owned_current_header")
    def _(f):
        owned(f)
        f.write(f.clinvar / "clinvar.hg38.vep.vcf.gz.installation-candidate.ready.json", "{}")
        return ["retire-owned", *common(f)]

    @case("retire_owned_explicit")
    def _(f):
        explicit = f.vcf(f.cache / "explicit/clinvar.custom.vcf.gz")
        f.save(clinvar_vcf=str(explicit))
        return ["retire-owned", *common(f)]

    @case("retire_owned_twice_with_source")
    def _(f):
        owned(f, '##VEP="v115" assembly="GRCh38"')
        pinned = f.vcf(f.cache / "pinned/raw.vcf.gz", header=None)
        f.vcf(f.root / "fresh.vcf.gz", header='##VEP="v114" assembly="GRCh38"')
        return [("retire-owned", *common(f), "--source-vcf", str(pinned)),
                ("copy", str(f.root / "fresh.vcf.gz"), str(f.target)),
                ["retire-owned", *common(f), "--source-vcf", str(pinned)]]

    def shards(f):
        names = ["chr1", "chr1_KI270706v1_random", "prefixchr1suffix"]
        for name in names:
            for ending in (".vcf.gz", ".vcf.gz.tbi", ".gnomAD.vcf.gz", ".vcf"):
                f.write(f.clinvar / f"clinvar.hg38.{name}{ending}", name + ending)
        return [value for name in names for value in ("--chromosome", name)]

    @case("retire_shards_after_retire_owned")
    def _(f):
        owned(f, '##VEP="v115" assembly="GRCh38"')
        chosen = shards(f)
        f.vcf(f.root / "fresh.vcf.gz")
        retire = ["retire-shards", *common(f), "--input-vcf", str(f.target), *chosen]
        return [("retire-owned", *common(f)),
                ("copy", str(f.root / "fresh.vcf.gz"), str(f.target)),
                tuple(retire), retire]

    @case("retire_shards_contigs_changed")
    def _(f):
        owned(f, '##VEP="v115" assembly="GRCh38"')
        chosen = shards(f)
        f.vcf(f.root / "fresh.vcf.gz")
        return [("retire-owned", *common(f)),
                ("copy", str(f.root / "fresh.vcf.gz"), str(f.target)),
                ("retire-shards", *common(f), "--input-vcf", str(f.target), *chosen),
                ["retire-shards", *common(f), "--input-vcf", str(f.target), "--chromosome", "chr2"]]

    @case("retire_shards_without_journal")
    def _(f):
        owned(f)
        return ["retire-shards", *common(f), "--input-vcf", str(f.target), *shards(f)]

    @case("retire_shards_other_input")
    def _(f):
        owned(f)
        return ["retire-shards", *common(f), "--input-vcf", str(f.clinvar / "other.vcf.gz"), "--chromosome", "chr1"]

    @case("retire_shards_unsafe_name")
    def _(f):
        owned(f, '##VEP="v115" assembly="GRCh38"')
        f.vcf(f.root / "fresh.vcf.gz")
        return [("retire-owned", *common(f)),
                ("copy", str(f.root / "fresh.vcf.gz"), str(f.target)),
                ["retire-shards", *common(f), "--input-vcf", str(f.target), "--chromosome", "a/b"]]

    # ---- install
    @case("install_explicit_preserved")
    def _(f):
        (f.source / "data").rename(f.root / "fixture_data")
        (f.source / "data").symlink_to(f.repo / "data")
        stats = f.cache / "stats/clinvar_2star_stats.pkl.gz"
        stats.parent.mkdir(parents=True)
        shutil.copyfile(f.repo / "data/ClinVar/clinvar_2star_stats.pkl.gz", stats)
        f.save(clinvar_gene_stat=str(stats))
        return ["install", *common(f)]

    @case("install_prerequisite_failure")
    def _(f):
        f.save()
        return ["install", *common(f)]

    # ---- argument errors
    @case("arguments_missing_config")
    def _(f):
        return ["plan", "--source-root", str(f.source)]

    @case("arguments_unknown_command")
    def _(f):
        return ["refresh", *common(f)]

    @case("arguments_retire_shards_without_chromosome")
    def _(f):
        return ["retire-shards", *common(f), "--input-vcf", str(f.target)]

    return cases


# ---------------------------------------------------------------- comparison

def normalize(text, archives=None):
    if not isinstance(text, str):
        return text
    for name, label in (archives or {}).items():
        text = text.replace(name, label)
    text = TIME.sub("<TIME>", text)
    text = ARCHIVE.sub("<ARCHIVE>", text)
    return RAW.sub(lambda m: f'"{m[1]}": "<ID>"', text)


def archive_labels(root):
    """Retained-generation folders named by their journal, not their random names."""
    labels = {}
    for folder, directories, _ in os.walk(root, followlinks=False):
        for name in directories:
            if ARCHIVE.fullmatch(name):
                journal = Path(folder) / name / "retirement.json"
                text = normalize(journal.read_text()) if journal.exists() else ""
                labels[name] = "<ARCHIVE:" + hashlib.sha256(text.encode()).hexdigest()[:12] + ">"
    return labels


def tree(root):
    archives = archive_labels(root)
    state = {}
    for folder, directories, names in os.walk(root, followlinks=False):
        for name in sorted(directories + names):
            path = Path(folder) / name
            relative = normalize(str(path.relative_to(root)), archives)
            if path.is_symlink():
                state[relative] = ["link", os.readlink(path)]
                continue
            mode = oct(path.lstat().st_mode & 0o7777)
            if path.is_dir():
                state[relative] = ["directory", mode]
                continue
            data = path.read_bytes()
            if len(data) <= 1 << 20:
                try:
                    state[relative] = ["text", normalize(data.decode("utf-8"), archives), mode]
                    continue
                except UnicodeDecodeError:
                    pass
            state[relative] = ["bytes", hashlib.sha256(data).hexdigest(), mode]
    return state


def run(command, folder, label, environment, cwd):
    with (folder / (label + ".stdout")).open("wb") as out, (folder / (label + ".stderr")).open("wb") as err:
        process = subprocess.run(list(map(str, command)), cwd=cwd, stdout=out, stderr=err,
                                 stdin=subprocess.DEVNULL, env=environment, timeout=600)
    return process.returncode


def outcome(folder, label, status, root):
    stdout = (folder / (label + ".stdout")).read_bytes().decode("utf-8", "replace")
    stderr = (folder / (label + ".stderr")).read_bytes().decode("utf-8", "replace")
    result = {"exit": status, "stdout": normalize(stdout), "tree": tree(root)}
    if label == "python":
        recorded = folder / "exception.json"
        value = json.loads(recorded.read_text()) if recorded.exists() else {}
        if value:
            result.update(category=value.get("category"), message=value.get("message"))
    else:
        matched = re.search(r"(?:\A|\n)priva: ([A-Za-z_][A-Za-z_0-9]*): ((?:(?!\npriva: ).)*)\n\Z", stderr, re.S)
        if matched:
            result.update(category=matched[1], message=matched[2])
        elif status not in (0, 2):
            result.update(category=None, message=stderr)
    if "message" in result:
        result["message"] = normalize(result["message"])
    return result


def compare(python, native, before):
    if python["exit"] == 2 or native["exit"] == 2:
        return python["exit"] == native["exit"] == 2 and python["tree"] == native["tree"] == before
    return python == native


def arm(args, producer, name, build, label, folder, python):
    """Build the case folder, run its setup steps, then the compared command."""
    root = (folder / "root").absolute()
    if root.exists():
        shutil.rmtree(root)
    root.mkdir()
    fixture = Fixture(root, args.repo)
    steps = build(fixture)
    fixture.finish()
    steps = steps if isinstance(steps[0], (list, tuple)) else [steps]
    environment = {key: value for key, value in os.environ.items() if not key.startswith(("PRIVA_", "FIXTURE_"))}
    environment.update(fixture.env, PYTHONDONTWRITEBYTECODE="1")
    cwd = args.work_root / "cwd"
    for step in steps[:-1]:
        if step[0] == "append":
            with open(step[1], "a") as handle:
                handle.write(step[2])
        elif step[0] == "copy":
            shutil.copyfile(step[1], step[2])
            os.utime(step[2], ns=(MTIME, MTIME))
        else:
            # Earlier steps of a case use the arm's own implementation.
            command = ([python, producer, *step] if label == "python"
                       else [args.binary, "clinvar-gene-refresh", *step])
            done = subprocess.run(list(map(str, command)), cwd=cwd, env=environment, capture_output=True, text=True)
            (folder / f"{label}.setup.log").open("a").write(done.stdout + done.stderr)
    before = tree(root)
    exception = folder / "exception.json"
    exception.unlink(missing_ok=True)
    if label == "python":
        command = [python, Path(__file__).resolve(), "oracle", "--producer", producer,
                   "--exception", exception, "--", *steps[-1]]
    else:
        command = [args.binary, "clinvar-gene-refresh", *steps[-1]]
    status = run(command, folder, label, environment, cwd)
    result = outcome(folder, label, status, root)
    root.rename(folder / f"{label}_root")
    return result, before


def check(args):
    args.repo = args.repo.absolute()
    producer = args.repo / "scripts/clinvar_gene_refresh.py"
    if hashlib.sha256(producer.read_bytes()).hexdigest() != PRODUCER_SHA256:
        raise SystemExit("clinvar_gene_refresh.py differs from the qualified source; update the pin")
    python = shutil.which("python")
    if not python:
        raise SystemExit("python is not on PATH")
    cases = catalogue()
    unknown = set(args.only or ()) - cases.keys()
    if unknown:
        raise SystemExit("Unknown cases: " + ", ".join(sorted(unknown)))
    args.work_root = args.work_root.absolute()
    args.work_root.mkdir(parents=True, exist_ok=False)
    (args.work_root / "cwd").mkdir()
    rows, skipped = {}, []
    for name, (build, network) in cases.items():
        if args.only and name not in args.only:
            continue
        if network and not args.network:
            skipped.append(name)
            continue
        folder = args.work_root / name
        folder.mkdir()
        results, befores = {}, {}
        for label in ("python", "native"):
            results[label], befores[label] = arm(args, producer, name, build, label, folder, python)
        passed = befores["python"] == befores["native"] and \
            compare(results["python"], results["native"], befores["native"])
        rows[name] = {"passed": passed, "python": results["python"], "native": results["native"]}
        write_json(folder / "comparison.json", rows[name])
        if passed:
            for label in ("python", "native"):
                shutil.rmtree(folder / f"{label}_root")
    failed = sorted(name for name, row in rows.items() if not row["passed"])
    report = {
        "producer": str(producer), "producer_sha256": PRODUCER_SHA256, "binary": str(args.binary),
        "binary_sha256": hashlib.sha256(Path(args.binary).read_bytes()).hexdigest(),
        "cases": len(rows), "python_successes": sum(row["python"]["exit"] == 0 for row in rows.values()),
        "skipped_network": skipped, "failed": failed, "rows": rows,
    }
    write_json(args.work_root / "complete.json", report)
    print(json.dumps({key: report[key] for key in ("cases", "python_successes", "skipped_network", "failed")}))
    if failed:
        raise SystemExit(1)


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    commands = parser.add_subparsers(dest="procedure", required=True)
    p = commands.add_parser("check")
    p.add_argument("--binary", required=True, type=Path)
    p.add_argument("--repo", required=True, type=Path)
    p.add_argument("--work-root", required=True, type=Path)
    p.add_argument("--network", action="store_true")
    p.add_argument("--only", nargs="*")
    p = commands.add_parser("oracle")
    p.add_argument("--producer", required=True, type=Path)
    p.add_argument("--exception", required=True, type=Path)
    p.add_argument("arguments", nargs=argparse.REMAINDER)
    args = parser.parse_args()
    if args.procedure == "oracle":
        if args.arguments[:1] == ["--"]:
            args.arguments = args.arguments[1:]
        oracle(args)
    else:
        check(args)


if __name__ == "__main__":
    main()
