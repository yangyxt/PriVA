#!/usr/bin/env python3
"""Compare `priva installer-prerequisites` with its unchanged Python source.

check --binary PRIVA --repo CHECKOUT --work-root NEW_DIRECTORY [--only CASE ...]
oracle --producer SCRIPT --exception JSON -- ARGUMENTS (internal reference launcher)

Each case builds two identical fixture folders. The Python arm runs
CHECKOUT/scripts/installer_prerequisites.py through the oracle launcher, which
records the class and text of the exception main() reports (or does not
catch) without changing main() or any check. The native arm runs
`PRIVA installer-prerequisites` with the same arguments. Compared exactly:
exit status, stdout (with the case folder written as <ARM>), the exception
category and message, and the complete fixture tree before and after (neither
arm may write anything, except that a case marked "writes" compares the
rewritten config.yaml text, with the case folder written as <ARM>, and its
permission bits, while every other entry must stay unchanged). A case with
"terminal" input runs both arms with standard input on a pseudo-terminal
without echo, holding that input; standard output stays a file, so both arms
write their prompts there. Declared differences, each named in its case:
- YAML syntax errors: PyYAML's subclass and text against the native YAMLError;
  only the category family is compared.
- decoder diagnostics inside "Cannot reuse bundled ClinVar gene statistics":
  the restricted pickle reader names other faults than pickle; the wrapper
  text and category are still compared.
- refused profiles: the native command must fail with the named category
  (anchors, --allow-pinned-reference-defaults, reducer-built pickles).
Writes WORK_ROOT/complete.json and exits 1 on any mismatch.
"""
from __future__ import annotations

import argparse
import gzip
import hashlib
import importlib.util
import json
import os
from pathlib import Path
import pickle
import pty
import re
import shutil
import subprocess
import sys
import termios

sys.dont_write_bytecode = True
PRODUCER_SHA256 = "3dc00df076de9a545e23e550696f85f4a77d2958c5c861ace9e9a1c955a86e6a"
YAML_ERRORS = {"YAMLError", "MarkedYAMLError", "ReaderError", "ScannerError", "ParserError",
               "ComposerError", "ConstructorError"}
CHR1 = {"hg19": 249250621, "hg38": 248956422}
MAPPING_HEADER = ("Gene stable ID\tTranscript stable ID\tHGNC symbol\t"
                  "UniProtKB/Swiss-Prot ID\tRefSeq mRNA ID\n")
MAPPING_ROW = "ENSG000001\tENST000001.4\tGENE\tP00001\tNM_000001.2\n"
GRCH37_MAPPING = "ensembl_GRCh37_mart116_nm_to_enst.tsv"
GRCH38_MAPPING = "ensembl_GRCh38_mart113_nm_to_enst.tsv"
FTP_MAPPING = "ensembl_GRCh38_mart116_ftp_candidate.tsv"
FIELDS = ("ensembl_id", "clinvar_pathogenic", "large_aachange", "small_aachange",
          "putative_nmd_variants")


def load(path, name):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def write_json(path, value):
    Path(path).write_text(json.dumps(value, sort_keys=True, indent=2) + "\n")


def oracle(args):
    """Run the unchanged main() and record the exception it reports."""
    module = load(args.producer, "original_installer_prerequisites")
    record = {}
    builtin_print = print

    def observed_print(*values, **options):
        # main() reports a caught failure with print(f"ERROR: {error}") inside
        # its except clause, where sys.exc_info() is that exception.
        if options.get("file") is sys.stderr and values and str(values[0]).startswith("ERROR: "):
            kind, value, _ = sys.exc_info()
            if kind is not None:
                record.update(category=kind.__name__, message=str(value))
        return builtin_print(*values, **options)

    module.print = observed_print
    sys.argv = [str(args.producer), *args.arguments]
    try:
        status = module.main()
    except SystemExit:
        raise
    except BaseException as exc:
        record.update(category=type(exc).__name__, message=str(exc), uncaught=True)
        write_json(args.exception, record)
        raise
    write_json(args.exception, record)
    raise SystemExit(status)


# ---------------------------------------------------------------- fixtures

def stats_pickle(protocol=4, value=None):
    if value is None:
        value = {"GENE": {field: (set() if field != "ensembl_id" else "ENSG1") for field in FIELDS}}
    return gzip.compress(pickle.dumps(value, protocol=protocol))


def archive(container=1, schema=5, extra=0):
    body = b"\0" * 16
    return (b"PRIVARKV" + container.to_bytes(4, "little") + schema.to_bytes(4, "little")
            + (1).to_bytes(4, "little") + (1).to_bytes(4, "little")
            + (len(body) + extra).to_bytes(8, "little") + body)


def chain(source=CHR1["hg19"], target=CHR1["hg38"]):
    return f"chain 1 chr1 {source} + 0 1 chr1 {target} + 0 1 1\n1\n"


def base_files():
    files = {
        "mapping/original_mapping.tsv": MAPPING_HEADER + MAPPING_ROW,
        "ref/hg19.fa": ">chr1\nA\n",
        "ref/hg19.fa.fai": f"chr1\t{CHR1['hg19']}\t6\t1\t2\n",
        "ref/hg38.fa": ">chr1\nA\n",
        "ref/hg38.fa.fai": f"chr1\t{CHR1['hg38']}\t6\t1\t2\n",
        "ref/hg19ToHg38.over.chain": chain(),
        "manual/primateai.tsv.gz": "fixture\n",
        "manual/spliceai_snv.prescore.gz": "fixture\n",
        "manual/spliceai_indel.prescore.gz": "fixture\n",
        "checkout/data/ClinVar/clinvar_2star_stats.pkl.gz": stats_pickle(),
        "checkout/data/alphamissense/alphamissense_mean_score.tsv": "fixture\n",
        "checkout/data/vep/csq_interpretations.tsv": "fixture\n",
        "checkout/data/hgnc/hgnc_non_alt_loci_set_sha256_890b18743cc8.tsv": "fixture\n",
        "checkout/data/dispensable_genes/dispensable_gene_list.txt": "fixture\n",
        "checkout/data/gofcards/gofcards_mechanism_reviews.tsv": "fixture\n",
        "checkout/data/mondo/disease_scope_overrides.tsv": "fixture\n",
        "checkout/data/patho_mechanism/gene_nonlof_mechanism_curated_assertions.schema.json": "{}\n",
        "cache/.keep": "",
        "tmp/.keep": "",
    }
    for assembly in ("hg19", "hg38"):
        files[f"checkout/data/HCSeeker/HC_spots.{assembly}.tsv"] = "fixture\n"
        files[f"checkout/data/gnomAD_RMC/gnomad_v2.1.1_rmc_flat.{assembly}.tsv"] = "fixture\n"
        files[f"checkout/data/repeats/ucsc.{assembly}.fasta.repeat.bed.gz"] = "fixture\n"
    for assembly in ("GRCh37", "GRCh38"):
        root = f"vep/homo_sapiens_merged/113_{assembly}/"
        files[root + "info.txt"] = f"species\thomo_sapiens\nassembly\t{assembly}\nvar_type\ttabix\n"
        for name in ("chr_synonyms.txt", "1/all_vars.gz", "1/all_vars.gz.csi"):
            files[root + name] = "fixture\n"
    return files


def base_config():
    return {
        "assembly": "hg38",
        "ref_genome": "{arm}/ref/hg38.fa",
        "nm_to_enst_tsv": "{arm}/mapping/original_mapping.tsv",
        "vep_installed_version": "113",
        "gofcards_hg19_fasta": "{arm}/ref/hg19.fa",
        "gofcards_hg38_fasta": "{arm}/ref/hg38.fa",
        "gofcards_hg19_to_hg38_chain": "{arm}/ref/hg19ToHg38.over.chain",
        "vep_cache_dir": "{arm}/vep",
        "cache_root": "{arm}/cache",
        "tmp_dir": "{arm}/tmp",
        "primateai_prescore": "{arm}/manual/primateai.tsv.gz",
        "spliceai_snv_prescore": "{arm}/manual/spliceai_snv.prescore.gz",
        "spliceai_indel_prescore": "{arm}/manual/spliceai_indel.prescore.gz",
    }


def yaml_value(value, arm):
    if value is None:
        return "null"
    if isinstance(value, bool):
        return "true" if value else "false"
    if isinstance(value, int):
        return str(value)
    if isinstance(value, list):
        return "[" + ", ".join(value) + "]"
    if isinstance(value, dict):
        return "".join(f"\n  {key}: {yaml_value(item, arm)}" for key, item in value.items())
    return json.dumps(value.replace("{arm}", str(arm)), ensure_ascii=False)


def config_text(case, arm):
    if "raw" in case:
        return case["raw"].replace("{arm}", str(arm))
    config = base_config()
    for key in case.get("unset", ()):
        config.pop(key, None)
    config.update(case.get("set", {}))
    lines = [f"{key}: {yaml_value(value, arm)}" for key, value in config.items()]
    return "# fixture configuration\n" + "\n".join(lines + case.get("append", [])) + "\n"


def materialize(arm, case):
    arm.mkdir(parents=True)
    files = base_files()
    files.update(case.get("files", {}))
    for relative, content in files.items():
        path = arm / relative
        if content is None:
            continue
        path.parent.mkdir(parents=True, exist_ok=True)
        if isinstance(content, tuple):
            path.symlink_to(content[1].replace("{arm}", str(arm)))
        elif isinstance(content, bytes):
            path.write_bytes(content)
        else:
            path.write_text(content.replace("{arm}", str(arm)))
    (arm / "config.yaml").write_bytes(case["config_bytes"] if "config_bytes" in case
                                      else config_text(case, arm).encode())
    if case.get("config_directory"):
        (arm / "config.yaml").unlink()
        (arm / "config.yaml").mkdir()


def tree(root):
    state = {}
    for folder, directories, names in os.walk(root, followlinks=False):
        for name in sorted(directories + names):
            path = Path(folder) / name
            relative = str(path.relative_to(root))
            if path.is_symlink():
                state[relative] = ["link", os.readlink(path)]
            elif path.is_dir():
                state[relative] = ["directory"]
            else:
                state[relative] = ["file", hashlib.sha256(path.read_bytes()).hexdigest()]
    return state


# ---------------------------------------------------------------- cases

def provenance(path, text):
    return {path: MAPPING_HEADER + MAPPING_ROW, path + ".provenance.json": text}


def catalogue():
    bundled38 = "checkout/data/hgnc/" + GRCH38_MAPPING
    cases = {
        # main stage
        "main_valid": {},
        "main_old_host": {"set": {"hpo_assertions": "/paedyl01/hpo.tsv.gz",
                                  "input_vcf": "/paedyl01/ignored.gz"}},
        "main_two_old_hosts": {"set": {"hpo_assertions": "/paedyl01/a", "clinvar_gene_stat": "/paedyl01/b"}},
        "main_mapping_missing": {"set": {"nm_to_enst_tsv": "{arm}/mapping/absent.tsv"}},
        "main_mapping_relative": {"set": {"nm_to_enst_tsv": "mapping/../mapping/local.tsv"},
                                  "files": {"checkout/mapping/local.tsv": MAPPING_HEADER + MAPPING_ROW}},
        "main_mapping_bad_header": {"files": {"mapping/original_mapping.tsv": "a\tb\tc\td\te\n" + MAPPING_ROW}},
        "main_mapping_short_header": {"files": {"mapping/original_mapping.tsv": "a\tb\n"}},
        "main_mapping_empty_header_line": {"files": {"mapping/original_mapping.tsv": "\n" + MAPPING_ROW}},
        "main_mapping_no_rows": {"files": {"mapping/original_mapping.tsv": MAPPING_HEADER + "x\ty\tz\tw\tv\n"}},
        "main_mapping_bom_crlf_quoted": {"files": {"mapping/original_mapping.tsv":
            "﻿" + MAPPING_HEADER.replace("\n", "\r\n").replace("RefSeq mRNA ID", " RefSeq mRNA ID ")
            + "\"a\tb\"\tc\n" + "G\t\"ENST2\"\tS\t\"multi\nline\"\tNM_2\r\n"}},
        "main_mapping_unterminated_quote": {"files": {"mapping/original_mapping.tsv":
            MAPPING_HEADER + "G\t\"ENST1\tS\tU\tNM_1\n"}},
        "main_mapping_invalid_utf8": {"files": {"mapping/original_mapping.tsv": b"Gene\xff\n"}},
        "main_mapping_field_limit": {"files": {"mapping/original_mapping.tsv": "x" * 131073 + "\n"}},
        "main_mapping_provenance_object": {"files": provenance("mapping/original_mapping.tsv", "{}\n")},
        "main_mapping_provenance_list": {"files": provenance("mapping/original_mapping.tsv", "[]\n")},
        "main_mapping_provenance_invalid_json": {"files": provenance("mapping/original_mapping.tsv", "{x\n")},
        "main_mapping_provenance_dangling": {"files": {"mapping/original_mapping.tsv.provenance.json":
                                                       ("link", "{arm}/absent.json")}},
        "main_mapping_pinned_name_no_sidecar": {
            "set": {"nm_to_enst_tsv": "{arm}/mapping/" + GRCH38_MAPPING},
            "files": {"mapping/" + GRCH38_MAPPING: MAPPING_HEADER + MAPPING_ROW}},
        "main_mapping_pinned_wrong_assembly": {
            "set": {"nm_to_enst_tsv": "{arm}/mapping/" + GRCH37_MAPPING},
            "files": {"mapping/" + GRCH37_MAPPING: MAPPING_HEADER + MAPPING_ROW}},
        "main_mapping_not_string": {"set": {"nm_to_enst_tsv": 12}},
        "main_default_bundled": {"set": {"nm_to_enst_tsv": None}, "files": {bundled38: "opaque bundle\n"}},
        "main_default_bundled_empty": {"set": {"nm_to_enst_tsv": ""}, "files": {bundled38: ""}},
        "main_default_bundled_link": {"set": {"nm_to_enst_tsv": None},
                                      "files": {bundled38: ("link", "{arm}/mapping/original_mapping.tsv")}},
        "main_default_missing_allowed": {"set": {"nm_to_enst_tsv": None},
                                         "args": ["--allow-missing-bundled-mapping"]},
        "main_default_missing_refused": {"set": {"nm_to_enst_tsv": "null"}},
        "main_default_orphan_provenance": {"set": {"nm_to_enst_tsv": None},
                                           "args": ["--allow-missing-bundled-mapping"],
                                           "files": {bundled38 + ".provenance.json": "{}\n"}},
        "main_default_ftp_116": {"set": {"nm_to_enst_tsv": None, "vep_installed_version": 116},
                                 "files": {"checkout/data/hgnc/" + FTP_MAPPING: "opaque bundle\n"}},
        "main_default_grch37": {"set": {"nm_to_enst_tsv": None, "assembly": "GRCh37"},
                                "files": {"checkout/data/hgnc/" + GRCH37_MAPPING: "opaque bundle\n"}},
        "main_vep_unsupported": {"set": {"nm_to_enst_tsv": None, "vep_installed_version": "999"}},
        "main_vep_flag_overrides": {"set": {"nm_to_enst_tsv": None}, "args": ["--vep-version", "999"]},
        "main_vep_bool": {"set": {"nm_to_enst_tsv": None, "vep_installed_version": True}},
        "main_vep_float": {"set": {"nm_to_enst_tsv": None}, "append": ["vep_installed_version: 113.0"]},
        "main_assembly_unknown": {"set": {"nm_to_enst_tsv": None, "assembly": "hg20"}},
        "main_assembly_missing": {"set": {"nm_to_enst_tsv": None}, "unset": ["assembly"]},
        "main_hpo_not_string": {"set": {"hpo_assertions": 5}},
        "main_hpo_blank": {"set": {"hpo_assertions": "   "}},
        "main_clinvar_missing_allowed": {"files": {"checkout/data/ClinVar/clinvar_2star_stats.pkl.gz": None},
                                         "args": ["--allow-missing-bundled-clinvar"]},
        "main_clinvar_empty_allowed": {"files": {"checkout/data/ClinVar/clinvar_2star_stats.pkl.gz": b""},
                                       "args": ["--allow-missing-bundled-clinvar"]},
        "main_clinvar_empty_refused": {"files": {"checkout/data/ClinVar/clinvar_2star_stats.pkl.gz": b""}},
        "main_clinvar_external_missing_allowed": {"set": {"clinvar_gene_stat": "{arm}/absent.pkl.gz"},
                                                  "args": ["--allow-missing-bundled-clinvar"]},
        "main_clinvar_link_not_deferred": {
            "files": {"checkout/data/ClinVar/clinvar_2star_stats.pkl.gz": ("link", "{arm}/absent.pkl.gz")},
            "args": ["--allow-missing-bundled-clinvar"]},
        "main_duplicate_keys_last_wins": {"append": ["assembly: hg19"],
                                          "set": {"nm_to_enst_tsv": None},
                                          "files": {"checkout/data/hgnc/" + GRCH37_MAPPING: "bundle\n"}},
        # ClinVar statistics stage
        "clinvar_bundled_valid": {"stage": "clinvar"},
        "clinvar_external_valid": {"stage": "clinvar", "set": {"clinvar_gene_stat": "{arm}/stats/x.pkl.gz"},
                                   "files": {"stats/x.pkl.gz": stats_pickle()}},
        "clinvar_external_protocol5": {"stage": "clinvar", "set": {"clinvar_gene_stat": "{arm}/stats/x.pkl.gz"},
                                       "files": {"stats/x.pkl.gz": stats_pickle(5)}},
        "clinvar_external_missing": {"stage": "clinvar", "set": {"clinvar_gene_stat": "{arm}/stats/x.pkl.gz"}},
        "clinvar_missing_fields": {"stage": "clinvar", "set": {"clinvar_gene_stat": "{arm}/stats/x.pkl.gz"},
                                   "files": {"stats/x.pkl.gz": stats_pickle(4, {"G": {"ensembl_id": "E"}})}},
        "clinvar_empty_dict": {"stage": "clinvar", "set": {"clinvar_gene_stat": "{arm}/stats/x.pkl.gz"},
                               "files": {"stats/x.pkl.gz": stats_pickle(4, {})}},
        "clinvar_not_dict": {"stage": "clinvar", "set": {"clinvar_gene_stat": "{arm}/stats/x.pkl.gz"},
                             "files": {"stats/x.pkl.gz": stats_pickle(4, ["a"])}},
        "clinvar_not_gzip": {"stage": "clinvar", "set": {"clinvar_gene_stat": "{arm}/stats/x.pkl.gz"},
                             "files": {"stats/x.pkl.gz": b"not gzip at all"}},
        "clinvar_truncated_gzip": {"stage": "clinvar", "loose_message": True,
                                   "set": {"clinvar_gene_stat": "{arm}/stats/x.pkl.gz"},
                                   "files": {"stats/x.pkl.gz": stats_pickle()[:20]}},
        "clinvar_global_pickle": {"stage": "clinvar", "loose_message": True,
                                  "set": {"clinvar_gene_stat": "{arm}/stats/x.pkl.gz"},
                                  "files": {"stats/x.pkl.gz": gzip.compress(pickle.dumps(
                                      __import__("collections").OrderedDict(a=1), protocol=4))}},
        "clinvar_protocol2_sets": {"stage": "clinvar", "refused": "ValueError",
                                   "set": {"clinvar_gene_stat": "{arm}/stats/x.pkl.gz"},
                                   "files": {"stats/x.pkl.gz": stats_pickle(2)}},
        "clinvar_rkyv_valid": {"stage": "clinvar", "set": {"clinvar_gene_stat": "{arm}/stats/x.rkyv"},
                               "files": {"stats/x.rkyv": archive()}},
        "clinvar_rkyv_schema": {"stage": "clinvar", "set": {"clinvar_gene_stat": "{arm}/stats/x.rkyv"},
                                "files": {"stats/x.rkyv": archive(schema=4)}},
        "clinvar_rkyv_length": {"stage": "clinvar", "set": {"clinvar_gene_stat": "{arm}/stats/x.rkyv"},
                                "files": {"stats/x.rkyv": archive(extra=1)}},
        "clinvar_rkyv_short": {"stage": "clinvar", "set": {"clinvar_gene_stat": "{arm}/stats/x.rkyv"},
                               "files": {"stats/x.rkyv": b"PRIVARKV"}},
        "clinvar_rkyv_magic": {"stage": "clinvar", "set": {"clinvar_gene_stat": "{arm}/stats/x.rkyv"},
                               "files": {"stats/x.rkyv": b"X" * 40}},
        "clinvar_relative_folder_missing": {"stage": "clinvar", "set": {"clinvar_gene_stat": "generated/r.pkl.gz"}},
        "clinvar_bare_name": {"stage": "clinvar", "set": {"clinvar_gene_stat": "release.pkl.gz"},
                              "files": {"checkout/data/ClinVar/release.pkl.gz": stats_pickle()}},
        "clinvar_integer_name": {"stage": "clinvar", "set": {"clinvar_gene_stat": 5}},
        "clinvar_false_is_default": {"stage": "clinvar", "set": {"clinvar_gene_stat": False}},
        "clinvar_bundled_empty": {"stage": "clinvar",
                                  "files": {"checkout/data/ClinVar/clinvar_2star_stats.pkl.gz": b""}},
        "clinvar_bundled_link_outside": {
            "stage": "clinvar",
            "files": {"stats/x.pkl.gz": stats_pickle(),
                      "checkout/data/ClinVar/clinvar_2star_stats.pkl.gz": ("link", "{arm}/stats/x.pkl.gz")}},
        # bundled snapshots
        "snapshots_valid": {"stage": "snapshots"},
        "snapshots_grch37": {"stage": "snapshots", "set": {"assembly": "GRCh37"}},
        "snapshots_bad_assembly": {"stage": "snapshots", "set": {"assembly": "hg20"}},
        "snapshots_list_assembly": {"stage": "snapshots", "set": {"assembly": ["hg38"]}},
        "snapshots_repeat_custom": {"stage": "snapshots", "set": {"repeat_region_file_name": "custom.bed.gz"},
                                    "files": {"checkout/data/repeats/custom.bed.gz": "fixture\n"}},
        "snapshots_repeat_path": {"stage": "snapshots", "set": {"repeat_region_file_name": "a/b.bed.gz"}},
        "snapshots_repeat_dot": {"stage": "snapshots", "set": {"repeat_region_file_name": ".."}},
        "snapshots_repeat_integer": {"stage": "snapshots", "set": {"repeat_region_file_name": 5}},
        "snapshots_missing_rmc": {"stage": "snapshots",
                                  "files": {"checkout/data/gnomAD_RMC/gnomad_v2.1.1_rmc_flat.hg38.tsv": None}},
        "snapshots_review_explicit_absent": {"stage": "snapshots",
                                             "set": {"gofcards_mechanism_review_tsv": "custom/review.tsv"}},
        "snapshots_review_explicit_default": {
            "stage": "snapshots",
            "set": {"gofcards_mechanism_review_tsv": "data/gofcards/gofcards_mechanism_reviews.tsv"},
            "files": {"checkout/data/gofcards/gofcards_mechanism_reviews.tsv": None}},
        "snapshots_hgnc_absolute_empty": {"stage": "snapshots", "set": {"hgnc_table": "{arm}/hgnc.tsv"},
                                          "files": {"hgnc.tsv": ""}},
        "snapshots_hgnc_absolute": {"stage": "snapshots", "set": {"hgnc_table": "{arm}/hgnc.tsv"},
                                    "files": {"hgnc.tsv": "x\n"}},
        # mapping stages
        "mapping_print": {"stage": "mapping", "args": ["--print-mapping"]},
        "mapping_quiet": {"stage": "mapping"},
        "mapping_target_default": {"stage": "mapping-target", "set": {"nm_to_enst_tsv": None}},
        "mapping_target_relative": {"stage": "mapping-target", "set": {"nm_to_enst_tsv": "./x//y.tsv"}},
        "mapping_target_not_string": {"stage": "mapping-target", "set": {"nm_to_enst_tsv": True}},
        "mapping_source_planned": {"stage": "mapping-source", "set": {"nm_to_enst_tsv": None}},
        "mapping_source_ftp": {"stage": "mapping-source", "set": {"nm_to_enst_tsv": None},
                               "args": ["--vep-version", "116"]},
        "mapping_source_explicit": {"stage": "mapping-source"},
        "mapping_source_dangling": {"stage": "mapping-source",
                                    "set": {"nm_to_enst_tsv": "{arm}/mapping/link.tsv"},
                                    "files": {"mapping/link.tsv": ("link", "{arm}/absent.tsv")}},
        "mapping_source_exists": {"stage": "mapping-source", "set": {"nm_to_enst_tsv": None},
                                  "files": {bundled38: "bundle\n"}},
        "mapping_source_orphan": {"stage": "mapping-source", "set": {"nm_to_enst_tsv": None},
                                  "files": {bundled38 + ".provenance.json": "{}\n"}},
        "mapping_workspace_cache_root": {"stage": "mapping-workspace"},
        "mapping_workspace_tmp": {"stage": "mapping-workspace", "set": {"cache_root": None}},
        "mapping_workspace_relative_inside": {"stage": "mapping-workspace", "set": {"cache_root": "rel/cache"}},
        "mapping_workspace_none": {"stage": "mapping-workspace", "set": {"cache_root": "", "tmp_dir": None}},
        "mapping_workspace_integer": {"stage": "mapping-workspace", "set": {"cache_root": 7}},
        # GoFCards stages
        "gofcards_valid": {"stage": "gofcards"},
        "gofcards_wrong_hg19": {"stage": "gofcards", "files": {"ref/hg19.fa.fai": "chr1\t1000\t6\t1\t2\n"}},
        "gofcards_same_fasta": {"stage": "gofcards", "set": {"gofcards_hg19_fasta": "{arm}/ref/hg38.fa"}},
        "gofcards_reverse_chain": {"stage": "gofcards",
                                   "files": {"ref/hg19ToHg38.over.chain": chain(CHR1["hg38"], CHR1["hg19"])}},
        "gofcards_chain_late_header": {"stage": "gofcards", "files": {"ref/hg19ToHg38.over.chain":
                                                                      "# comment\r\nchainx\n" + chain()}},
        "gofcards_chain_no_header": {"stage": "gofcards", "files": {"ref/hg19ToHg38.over.chain": "1\n"}},
        "gofcards_chain_noninteger": {"stage": "gofcards", "files": {"ref/hg19ToHg38.over.chain":
                                                                     chain().replace("249250621", "x")}},
        "gofcards_chain_gzip": {"stage": "gofcards",
                                "set": {"gofcards_hg19_to_hg38_chain": "{arm}/ref/forward.chain.gz"},
                                "files": {"ref/forward.chain.gz": gzip.compress(chain().encode())}},
        "gofcards_chain_bad_gzip": {"stage": "gofcards",
                                    "set": {"gofcards_hg19_to_hg38_chain": "{arm}/ref/forward.chain.gz"},
                                    "files": {"ref/forward.chain.gz": chain()}},
        "gofcards_chain_missing": {"stage": "gofcards", "unset": ["gofcards_hg19_to_hg38_chain"]},
        "gofcards_fai_unicode_digits": {"stage": "gofcards", "files": {"ref/hg38.fa.fai":
                                                                       "chr1\t２４８９５６４２２\t6\t1\t2\n"}},
        "gofcards_fai_bad_integer": {"stage": "gofcards", "files": {"ref/hg38.fa.fai": "chr2\tlong\n"}},
        "gofcards_fai_alias_and_duplicates": {"stage": "gofcards", "files": {
            "ref/hg38.fa.fai": f"1\t{CHR1['hg38']}\nchr1\t5\nchr1\t{CHR1['hg38']} \nshort\n"}},
        "gofcards_fai_missing": {"stage": "gofcards", "files": {"ref/hg19.fa.fai": None}},
        "gofcards_cache_missing": {"stage": "gofcards",
                                   "files": {"vep/homo_sapiens_merged/113_GRCh38/info.txt": None}},
        "gofcards_cache_no_index": {"stage": "gofcards",
                                    "files": {"vep/homo_sapiens_merged/113_GRCh37/1/all_vars.gz.csi": None}},
        "gofcards_cache_tbi_only": {"stage": "gofcards", "files": {
            "vep/homo_sapiens_merged/113_GRCh37/1/all_vars.gz.csi": None,
            "vep/homo_sapiens_merged/113_GRCh37/1/all_vars.gz.tbi": "fixture\n"}},
        "gofcards_cache_wrong_metadata": {"stage": "gofcards", "files": {
            "vep/homo_sapiens_merged/113_GRCh38/info.txt": "species homo_sapiens\nassembly GRCh37\nvar_type tabix\n"}},
        "gofcards_cache_info_spacing": {"stage": "gofcards", "files": {
            "vep/homo_sapiens_merged/113_GRCh38/info.txt":
                "  species   homo_sapiens  \r\nassembly\tGRCh38\x0bvar_type tabix\nlonely\n"}},
        "gofcards_version_flag": {"stage": "gofcards", "args": ["--vep-version", "116"]},
        "gofcards_version_nondigit": {"stage": "gofcards", "set": {"vep_installed_version": "11a"}},
        "gofcards_version_missing": {"stage": "gofcards", "unset": ["vep_installed_version"]},
        "gofcards_caches_deferred": {"stage": "gofcards", "args": ["--allow-missing-merged-caches",
                                                                   "--vep-version", "116"]},
        "gofcards_caches_deferred_smoke": {"stage": "gofcards", "env": {"PRIVA_SMOKE_SKIP_VEP": "1"},
                                           "args": ["--allow-missing-merged-caches", "--vep-version", "116"]},
        "gofcards_cache_parent_file": {"stage": "gofcards", "args": ["--allow-missing-merged-caches"],
                                       "set": {"vep_cache_dir": "{arm}/vepfile"}, "files": {"vepfile": "x\n"}},
        "gofcards_main_selected_reference_pending": {
            "stage": "gofcards-main", "args": ["--allow-missing-merged-caches", "--vep-version", "116"],
            "set": {"cache_root": None}, "files": {"ref/hg38.fa": None, "ref/hg38.fa.fai": None}},
        "gofcards_main_wrong_existing_reference": {
            "stage": "gofcards-main", "set": {"ref_genome": "{arm}/ref/hg19.fa",
                                              "gofcards_hg38_fasta": "{arm}/ref/hg19.fa"}},
        "gofcards_main_cache_root": {"stage": "gofcards-main", "set": {"cache_folders": {"vep_cache": "myvep"}}},
        "gofcards_main_cache_root_default": {"stage": "gofcards-main", "files": {
            "cache/VEP_cache/homo_sapiens_merged/113_GRCh37/info.txt": None}},
        "gofcards_main_cache_folders_list": {"stage": "gofcards-main", "set": {"cache_folders": ["a"]}},
        "gofcards_main_cache_folders_null_value": {"stage": "gofcards-main",
                                                   "set": {"cache_folders": {"vep_cache": None}}},
        "gofcards_main_cache_root_integer": {"stage": "gofcards-main", "set": {"cache_root": 3}},
        "gofcards_pinned_defaults": {"stage": "gofcards", "refused": "NotImplementedError",
                                     "args": ["--allow-pinned-reference-defaults"]},
        "gofcards_list_assembly": {"stage": "gofcards", "set": {"assembly": ["x"]}},
        # manual publisher inputs
        "manual_prepared": {"stage": "manual-inputs"},
        "manual_prepared_prompt": {"stage": "manual-inputs", "args": ["--prompt"]},
        "manual_missing": {"stage": "manual-inputs", "set": {"primateai_prescore": None,
                                                             "spliceai_indel_prescore": " null "}},
        "manual_missing_prompt_without_terminal": {"stage": "manual-inputs", "args": ["--prompt"],
                                                   "set": {"spliceai_snv_prescore": ""}},
        "manual_staging_name": {"stage": "manual-inputs",
                                "set": {"primateai_prescore": "{arm}/manual/primateai.tsv.gz.partial"},
                                "files": {"manual/primateai.tsv.gz.partial": "x\n"}},
        "manual_link_to_staging": {"stage": "manual-inputs", "files": {
            "manual/.download.primateai": "x\n",
            "manual/primateai.tsv.gz": ("link", "{arm}/manual/.download.primateai")}},
        "manual_empty_file": {"stage": "manual-inputs", "files": {"manual/primateai.tsv.gz": ""}},
        "manual_integer": {"stage": "manual-inputs", "set": {"primateai_prescore": 4}},
        "manual_duplicate_key": {"stage": "manual-inputs", "append": ["assembly: hg19"]},
        "manual_nested_duplicate_key": {"stage": "manual-inputs", "append": ["nested:", "  a: 1", "  a: 2"]},
        "manual_not_mapping": {"stage": "manual-inputs", "raw": "- a\n- b\n"},
        "manual_missing_config": {"stage": "manual-inputs", "config_arg": "absent/config.yaml"},
        # manual publisher inputs collected at a terminal
        "manual_prompt_fills_unset": {
            "stage": "manual-inputs", "args": ["--prompt"], "writes": True,
            "set": {"primateai_prescore": None, "spliceai_indel_prescore": " null "},
            "terminal": "{arm}/manual/primateai.tsv.gz\n  {arm}/manual/spliceai_indel.prescore.gz \t\n"},
        "manual_prompt_appends_absent_key": {
            "stage": "manual-inputs", "args": ["--prompt"], "writes": True,
            "unset": ["spliceai_snv_prescore"], "terminal": "{arm}/manual/spliceai_snv.prescore.gz\n"},
        "manual_prompt_keeps_comment": {
            "stage": "manual-inputs", "args": ["--prompt"], "writes": True, "unset": ["primateai_prescore"],
            "append": ["primateai_prescore: null   # fill at the terminal"],
            "terminal": "{arm}/manual/primateai.tsv.gz\n"},
        "manual_prompt_unicode_path": {
            "stage": "manual-inputs", "args": ["--prompt"], "writes": True,
            "set": {"primateai_prescore": ""}, "files": {"manual/primateai_é.tsv.gz": "fixture\n"},
            "terminal": "{arm}/manual/primateai_é.tsv.gz\n"},
        "manual_prompt_invalid_answer": {
            "stage": "manual-inputs", "args": ["--prompt"], "set": {"primateai_prescore": None},
            "terminal": "{arm}/manual/absent.tsv.gz\n"},
        "manual_prompt_later_answer_invalid": {
            "stage": "manual-inputs", "args": ["--prompt"],
            "set": {"primateai_prescore": None, "spliceai_snv_prescore": None},
            "terminal": "{arm}/manual/primateai.tsv.gz\n{arm}/manual/primateai.tsv.gz.partial\n",
            "files": {"manual/primateai.tsv.gz.partial": "x\n"}},
        "manual_prompt_end_of_input": {
            "stage": "manual-inputs", "args": ["--prompt"],
            "set": {"primateai_prescore": None, "spliceai_snv_prescore": None},
            "terminal": "{arm}/manual/primateai.tsv.gz\n\x04"},
        "manual_prompt_quoted_key": {
            "stage": "manual-inputs", "args": ["--prompt"], "unset": ["primateai_prescore"],
            "append": ['"primateai_prescore": null'], "terminal": "{arm}/manual/primateai.tsv.gz\n"},
        # restrictions, admission and configuration reading
        "restrict_prompt": {"args": ["--prompt"]},
        "restrict_clinvar_flag": {"stage": "snapshots", "args": ["--allow-missing-bundled-clinvar"]},
        "restrict_pinned_flag": {"stage": "mapping", "args": ["--allow-pinned-reference-defaults"]},
        "restrict_merged_flag": {"args": ["--allow-missing-merged-caches"]},
        "admission_invalid_stage": {"stage": "everything"},
        "admission_explicit_flag_value": {"args": ["--prompt=yes"]},
        "config_missing": {"config_arg": "absent.yaml"},
        "config_directory": {"config_directory": True},
        "config_empty": {"raw": "# nothing\n"},
        "config_braces": {"raw": "{}\n"},
        "config_list": {"raw": "[]\n"},
        "config_invalid_utf8": {"config_bytes": b"assembly: hg38\n\xff\n"},
        "config_tab": {"raw": "assembly:\thg38\n\tx: 1\n", "yaml_error": True},
        "config_anchor": {"raw": "assembly: &a hg38\n", "refused": "NotImplementedError"},
        "config_relative_argument": {"config_arg": "config.yaml"},
        "config_crlf": {"crlf": True},
        # inspection-root policy
        "policy_inside": {"env": {"PRIVA_INSTALL_INSPECT_ROOTS": "{arm}"}},
        "policy_outside": {"env": {"PRIVA_INSTALL_INSPECT_ROOTS": "{arm}/checkout:{arm}/ref"}},
        "policy_link_escape": {"env": {"PRIVA_INSTALL_INSPECT_ROOTS": "{arm}/mapping:{arm}/checkout"},
                               "set": {"nm_to_enst_tsv": "{arm}/mapping/linked.tsv"},
                               "files": {"outside/real.tsv": MAPPING_HEADER + MAPPING_ROW,
                                         "mapping/linked.tsv": ("link", "{arm}/outside/real.tsv")}},
        "policy_link_between_roots": {"env": {"PRIVA_INSTALL_INSPECT_ROOTS": "{arm}/mapping:{arm}/outside:{arm}/checkout"},
                                      "set": {"nm_to_enst_tsv": "{arm}/mapping/linked.tsv"},
                                      "files": {"outside/real.tsv": MAPPING_HEADER + MAPPING_ROW,
                                                "mapping/linked.tsv": ("link", "../outside/real.tsv")}},
        "policy_nested_roots": {"stage": "gofcards", "env": {"PRIVA_INSTALL_INSPECT_ROOTS": "{arm}:{arm}/vep"}},
    }
    for name, case in cases.items():
        case.setdefault("stage", "main")
    return cases


# ---------------------------------------------------------------- runs

def arguments(case, arm):
    config = case.get("config_arg", str(arm / "config.yaml"))
    return [config, "--source-root", str(arm / "checkout"), "--stage", case["stage"],
            *[value.replace("{arm}", str(arm)) for value in case.get("args", [])]]


def run(command, folder, label, environment, terminal=None):
    with (folder / (label + ".stdout")).open("wb") as out, (folder / (label + ".stderr")).open("wb") as err:
        if terminal is None:
            process = subprocess.run(list(map(str, command)), cwd=folder / label, stdout=out, stderr=err,
                                     stdin=subprocess.DEVNULL, env=environment, timeout=120)
            return process.returncode
        # Canonical mode keeps line reads and turns a leading Ctrl-D into end of input.
        controller, terminal_end = pty.openpty()
        try:
            attributes = termios.tcgetattr(terminal_end)
            attributes[3] &= ~termios.ECHO
            termios.tcsetattr(terminal_end, termios.TCSANOW, attributes)
            os.write(controller, terminal)
            process = subprocess.run(list(map(str, command)), cwd=folder / label, stdout=out, stderr=err,
                                     stdin=terminal_end, env=environment, timeout=120)
        finally:
            os.close(terminal_end)
            os.close(controller)
    return process.returncode


def normalize(text, arm):
    if not isinstance(text, str):
        return text
    for spelling in sorted({str(arm), os.path.realpath(arm)}, key=len, reverse=True):
        text = text.replace(spelling, "<ARM>")
    return text


def outcome(folder, label, status, arm, before, case):
    stdout = (folder / (label + ".stdout")).read_bytes().decode("utf-8", "replace")
    stderr = (folder / (label + ".stderr")).read_bytes().decode("utf-8", "replace")
    after = tree(arm)
    result = {"exit": status, "stdout": normalize(stdout, arm), "unchanged": after == before}
    if case.get("writes"):
        config = arm / "config.yaml"
        others = [{key: value for key, value in state.items() if key != "config.yaml"}
                  for state in (before, after)]
        result.update(others_unchanged=others[0] == others[1],
                      config=normalize(config.read_bytes().decode("utf-8", "replace"), arm),
                      config_mode=oct(config.stat().st_mode & 0o7777))
    if label == "python":
        recorded = folder / "exception.json"
        # A successful run records an empty object: no exception to compare.
        value = json.loads(recorded.read_text()) if recorded.exists() else {}
        if value:
            result.update(category=value.get("category"), message=value.get("message"))
    else:
        matched = re.fullmatch(r"priva: ([A-Za-z_][A-Za-z_0-9]*): (.*)\n", stderr, re.S)
        if matched:
            result.update(category=matched[1], message=matched[2])
        elif status not in (0, 2):
            result.update(category=None, message=stderr)
    if "message" in result:
        result["message"] = normalize(result["message"], arm)
    return result


def compare(case, python, native):
    if case.get("refused"):
        return native["exit"] == 1 and native.get("category") == case["refused"] and native["unchanged"]
    if python["exit"] == 2 or native["exit"] == 2:
        # Argument admission: argparse usage text differs by design.
        return python["exit"] == native["exit"] == 2 and python["unchanged"] and native["unchanged"]
    left, right = dict(python), dict(native)
    if left.get("category") in YAML_ERRORS and right.get("category") == "YAMLError":
        left["category"] = "YAMLError"
        left.pop("message"), right.pop("message")
    if case.get("loose_message") and "message" in left and "message" in right:
        prefix = re.compile(r"(Cannot reuse bundled ClinVar gene statistics [^:]*): .*?(\. Restore a .*)", re.S)
        left["message"], right["message"] = (prefix.sub(r"\1: <DECODER>\2", left["message"]),
                                             prefix.sub(r"\1: <DECODER>\2", right["message"]))
    left.pop("uncaught", None)
    if case.get("writes"):
        # Both arms rewrite config.yaml; the text and mode are compared above.
        return left == right and left["others_unchanged"]
    return left == right and left["unchanged"]


def check(args):
    producer = args.repo / "scripts/installer_prerequisites.py"
    if hashlib.sha256(producer.read_bytes()).hexdigest() != PRODUCER_SHA256:
        raise SystemExit("installer_prerequisites.py differs from the qualified source; update the pin")
    cases = catalogue()
    unknown = set(args.only or ()) - cases.keys()
    if unknown:
        raise SystemExit("Unknown cases: " + ", ".join(sorted(unknown)))
    args.work_root.mkdir(parents=True, exist_ok=False)
    rows = {}
    for name, case in cases.items():
        if args.only and name not in args.only:
            continue
        folder = args.work_root / name
        folder.mkdir()
        results = {}
        for label in ("python", "native"):
            arm = (folder / label).absolute()
            materialize(arm, case)
            if case.get("crlf"):
                path = arm / "config.yaml"
                path.write_bytes(path.read_bytes().replace(b"\n", b"\r\n"))
            before = tree(arm)
            environment = {key: value for key, value in os.environ.items()
                           if key not in ("PRIVA_INSTALL_INSPECT_ROOTS", "PRIVA_SMOKE_SKIP_VEP")}
            environment.update({key: value.replace("{arm}", str(arm))
                                for key, value in case.get("env", {}).items()})
            environment["PYTHONDONTWRITEBYTECODE"] = "1"
            if label == "python":
                command = [sys.executable, Path(__file__).resolve(), "oracle", "--producer", producer,
                           "--exception", folder / "exception.json", "--", *arguments(case, arm)]
            else:
                command = [args.binary, "installer-prerequisites", *arguments(case, arm)]
            terminal = case.get("terminal")
            if terminal is not None:
                terminal = terminal.replace("{arm}", str(arm)).encode()
            status = run(command, folder, label, environment, terminal)
            results[label] = outcome(folder, label, status, arm, before, case)
        passed = compare(case, results["python"], results["native"])
        rows[name] = {"passed": passed, "python": results["python"], "native": results["native"],
                      "refused": case.get("refused"), "loose_message": bool(case.get("loose_message"))}
        write_json(folder / "comparison.json", rows[name])
    failed = sorted(name for name, row in rows.items() if not row["passed"])
    report = {
        "producer": str(producer), "producer_sha256": PRODUCER_SHA256,
        "binary": str(args.binary),
        "binary_sha256": hashlib.sha256(Path(args.binary).read_bytes()).hexdigest(),
        "cases": len(rows), "python_successes": sum(row["python"]["exit"] == 0 for row in rows.values()),
        "failed": failed, "rows": rows,
    }
    write_json(args.work_root / "complete.json", report)
    print(json.dumps({key: report[key] for key in ("cases", "python_successes", "failed")}))
    if failed:
        raise SystemExit(1)


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    commands = parser.add_subparsers(dest="procedure", required=True)
    p = commands.add_parser("check")
    p.add_argument("--binary", type=Path, required=True)
    p.add_argument("--repo", type=Path, required=True)
    p.add_argument("--work-root", type=Path, required=True)
    p.add_argument("--only", action="append")
    p.set_defaults(function=check)
    p = commands.add_parser("oracle")
    p.add_argument("--producer", type=Path, required=True)
    p.add_argument("--exception", type=Path, required=True)
    p.add_argument("arguments", nargs=argparse.REMAINDER)
    p.set_defaults(function=oracle)
    args = parser.parse_args()
    if args.procedure == "oracle" and args.arguments[:1] == ["--"]:
        args.arguments = args.arguments[1:]
    args.function(args)


if __name__ == "__main__":
    main()
