#!/usr/bin/env python3
"""Generate a separate installation config without touching the source config."""
import argparse
from pathlib import Path
import re
import yaml


def transform(config, source_root, cache_root, data_root, threads):
    config = dict(config)
    vep_release = str(config.get("vep_installed_version", ""))
    if not re.fullmatch(r"[1-9][0-9]*", vep_release):
        raise ValueError("vep_installed_version must be a numeric VEP major release")
    for key, value in config.items():
        if isinstance(value, str) and value.startswith("/paedyl01/disk1/yangyxt/PriVA/"):
            config[key] = str(source_root / value.split("/PriVA/", 1)[1])
    config.update({
        "assembly": "hg19", "tmp_dir": str(data_root / "tmp"),
        "ref_genome": str(cache_root / "references/ucsc.hg19.fasta"),
        "cache_root": str(cache_root), "vep_installed_version": vep_release,
        "input_vcf": "", "ped_file": "", "alt_disease_vcf": "",
        "control_vcf": "", "pp1_vcf": "", "pp1_ped": "",
        "output_dir": str(data_root / "results"),
        "hub_file_dir": str(data_root / "annotation_hub/hg19"),
        "threads": threads, "threads_per_fam": min(threads, 4),
        "gene_mechanism_raw_dir": str(cache_root / "PriVA_mechanism_raw"),
        "clinvar_vcv_xml": str(cache_root / "PriVA_mechanism_raw/clinvar_vcv/ClinVarVCVRelease_00-latest_weekly.xml.gz"),
        "gofcards_hg38_fasta": str(cache_root / "references/ucsc.hg38.fasta"),
        "gofcards_hg19_to_hg38_chain": str(cache_root / "references/hg19ToHg38.over.chain"),
    })
    folders = {"vep_cache_dir": "VEP_cache", "vep_plugins_dir": "VEP_plugin",
               "vep_plugins_cachedir": "VEP_plugin_cache", "gnomad_vcf_dir": "gnomAD",
               "clinvar_vcf_dir": "ClinVar", "cadd_parent_dir": "CADD"}
    for key, folder in folders.items():
        old = config.get(key, "")
        new = str(cache_root / folder)
        if old:
            for child, value in list(config.items()):
                if isinstance(value, str) and value.startswith(old + "/"):
                    config[child] = new + value[len(old):]
        config[key] = new
    for key, value in config.items():
        if isinstance(value, str) and value.startswith("/paedyl01/"):
            raise ValueError(f"Unresolved old-host path {key}: {value}")
    return config


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--template", type=Path, required=True)
    parser.add_argument("--source-root", type=Path, required=True)
    parser.add_argument("--cache-root", type=Path, required=True)
    parser.add_argument("--data-root", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--threads", type=int, default=4)
    args = parser.parse_args()
    if args.output.exists():
        parser.error("Refusing to overwrite an existing installation config")
    result = transform(yaml.safe_load(args.template.read_text()), args.source_root,
                       args.cache_root, args.data_root, args.threads)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text("# Cluster installation profile; family inputs will be supplied later.\n"
                           "# Original scientific settings and source config remain preserved.\n"
                           + yaml.safe_dump(result, sort_keys=False))
    print(args.output)


if __name__ == "__main__":
    main()
