#!/usr/bin/env python3
"""Report installation gaps without treating downloads as scientific qualification."""
import argparse
from datetime import datetime, timezone
import json
from pathlib import Path
import re
import yaml

RUNTIME_KEYS = (
    "ref_genome", "gnomad_vcf_chrX", "clinvar_vcf", "clinvar_aa_stat",
    "clinvar_splice_stat", "clinvar_patho_af_stat", "clinvar_patho_exon_af_stat",
    "utr_annotator_file", "human_ancestor_fasta", "loftee_conservation_file",
    "alphamissense_prescore", "alphamissense_vep_vcf", "alphamissense_pd_stat",
    "alphamissense_tranx_domain_map", "alphamissense_intolerant_domains",
    "alphamissense_intolerant_motifs", "all_intolerant_domains", "pm1_regions_pkl",
    "conservation_file", "splicevault_prescore", "primateai_prescore",
    "spliceai_snv_prescore", "spliceai_indel_prescore", "cds_fasta_file",
    "clingen_map", "gene_dosage_sensitivity", "gene_nonlof_mechanism_json",
    "hpo_condition_mechanism_json", "hpo_assertions", "gofcards_exact_gof_cache",
    "mondo_disease_scope_registry", "nm_to_enst_tsv",
)


def inspect(path):
    if not path:
        return "not_configured"
    path = Path(path)
    if not path.is_file():
        return "missing"
    if path.stat().st_size == 0:
        return "empty"
    with path.open("rb") as stream:
        if stream.read(64).startswith(b"version https://git-lfs.github.com/spec/"):
            return "git_lfs_pointer"
    return "present_not_scientifically_qualified"


def licensed_blockers(resources, assembly="GRCh37"):
    """Report each missing configured cache independently; presence is not qualification."""
    statuses = {item["key"]: item["status"] for item in resources}
    splice_label = "SpliceAI2" if assembly == "GRCh38" else "SpliceAI v1"
    splice_action = ("Register the completed GRCh38 SpliceAI2 VCF/index pair and adjacent adapter receipt through the installer. Do not use partials or fall back to v1."
                     if assembly == "GRCh38" else
                     "Register the confirmed matching-assembly SNV/indel originals from caches/SpliceAI_v1 through the real installer after integrity checks. Do not duplicate the browser acquisition.")
    groups = (
        (splice_label, ("spliceai_snv_prescore", "spliceai_indel_prescore"),
         splice_action),
        ("PrimateAI v0.2", ("primateai_prescore",),
         "Both licensed originals were acquired. Import the confirmed matching-assembly file from caches/PrimateAI_v0_2 through the installer; validate its formatted/indexed view. Do not duplicate acquisition. A missing whole cache is not an ordinary missing variant score."),
    )
    return [{"resource": label, "action": action, "unavailable_config_keys": missing}
            for label, keys, action in groups
            if (missing := [key for key in keys
                            if statuses.get(key) != "present_not_scientifically_qualified"])]


def environment_receipt_matches(receipt, requested_release):
    """An earlier VEP release's software receipt cannot qualify this profile."""
    if not isinstance(receipt, dict) or receipt.get("status") != "software_checked":
        return False
    packages = receipt.get("packages")
    if not isinstance(packages, dict):
        return False
    version = re.fullmatch(r"([0-9]+)(?:\.[0-9]+)*", str(packages.get("ensembl-vep", "")))
    return version is not None and version.group(1) == str(requested_release)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", type=Path, required=True)
    parser.add_argument("--state-dir", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    config = yaml.safe_load(args.config.read_text())
    vep_release = str(config.get("vep_installed_version", ""))
    if not re.fullmatch(r"[1-9][0-9]*", vep_release):
        raise ValueError("vep_installed_version must be a numeric VEP major release")
    assembly = {"hg19": "GRCh37", "GRCh37": "GRCh37",
                "hg38": "GRCh38", "GRCh38": "GRCh38"}.get(config.get("assembly"))
    if assembly is None:
        raise ValueError("assembly must be hg19/GRCh37 or hg38/GRCh38")
    resources = [{"key": key, "path": config.get(key), "status": inspect(config.get(key))} for key in RUNTIME_KEYS]
    source_root = Path(config["gene_mechanism_builder_script"]).parents[1]
    selected_gene_stat = str(config.get("clinvar_gene_stat") or "")
    gene_stat_path = (Path(selected_gene_stat).absolute() if "/" in selected_gene_stat
                      else source_root / "data/ClinVar" / selected_gene_stat)
    relative_resources = {
        "clinvar_gene_stat_used_by_shell": gene_stat_path,
        "interpro_mapping_pickle": source_root / "data/InterPro" / Path(config["interpro_mapping_pickle"]).name,
        "loeuf_used_by_shell": source_root / "data/loeuf" / Path(config["loeuf_prescore"]).name,
        "repeat_regions_used_by_shell": source_root / "data/repeats" / Path(config["repeat_region_file_name"]).name,
        "alphamissense_mean_used_by_shell": source_root / "data/alphamissense/alphamissense_mean_score.tsv",
        "hgnc_used_by_scientific_code": source_root / "data/hgnc/non_alt_loci_set.tsv",
        "dispensable_genes_default": source_root / "data/dispensable_genes/dispensable_gene_list.txt",
    }
    for key, path in relative_resources.items():
        resources.append({"key": key, "path": str(path), "status": inspect(path)})
    if assembly == "GRCh38":
        resources.append({"key": "mavedb_file", "path": config.get("mavedb_file"), "status": inspect(config.get("mavedb_file"))})
    if config.get("pext_enabled"):
        resources.append({"key": "pext_mean_bw", "path": config.get("pext_mean_bw"), "status": inspect(config.get("pext_mean_bw"))})
    vep_info = Path(config["vep_cache_dir"]) / "homo_sapiens_merged" / f"{vep_release}_{assembly}" / "info.txt"
    resources.append({"key": "vep_cache_info", "path": str(vep_info), "status": inspect(vep_info)})
    receipts = {}
    for path in sorted(args.state_dir.glob("*.ready.json")):
        receipts[path.name] = json.loads(path.read_text())
    report = {
        "schema_version": 1, "checked_utc": datetime.now(timezone.utc).isoformat(),
        "status": "incomplete_installation", "config": str(args.config),
        "software_qualified": environment_receipt_matches(receipts.get("environment.ready.json"), vep_release),
        "resources": resources, "completed_stage_receipts": receipts,
        "approved_omissions": ["MaveDB GRCh37 (explicit user approval)"] if assembly == "GRCh37" and not config.get("mavedb_file") else [],
        "blockers": licensed_blockers(resources, assembly) + (
            [{"resource": "MaveDB hg19 configuration", "action": "The approved GRCh37 omission requires an empty mavedb_file; do not use a GRCh38 cache."}]
            if assembly == "GRCh37" and config.get("mavedb_file") else []
        ) + [
            {"resource": f"Shared PriVA/CADD official VEP{vep_release}", "action": "Qualify one annotation pass, transcript/peptide/GTF views and complete feature mapping. Packaged CADD VEP110 is reference-only."},
            {"resource": "Derived calibration caches", "action": "Complete matching-assembly population resources and AlphaMissense/ClinVar/domain/PM1 builders. Pass the acquired pinned BioMart mapping explicitly and validate coverage against the selected transcript release. It is a new candidate, not identical to the unavailable original mapping; do not substitute MANE. GRCh37 gnomAD liftover remains separate from GRCh38 source deployment."},
            {"resource": "Family acceptance", "action": "BH14233 GRCh38 packet is delivered and intake checked. Validate nonjoint multi-sample import, controls and unknown sex/affection handling, then compare exact variant ordering and ACMG judgments within one completed scientific profile. Public-data acceptance also remains pending."},
        ],
        "note": "Presence checks and completed download receipts are not end-to-end validation. Run this command on the host, outside filename-filtered AI views.",
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    temporary = args.output.with_name(args.output.name + ".partial")
    temporary.write_text(json.dumps(report, indent=2) + "\n")
    temporary.replace(args.output)
    missing = sum(item["status"] != "present_not_scientifically_qualified" for item in resources)
    print(f"Software qualified: {report['software_qualified']}; {missing}/{len(resources)} configured resource paths unavailable; installation incomplete.")


if __name__ == "__main__":
    main()
