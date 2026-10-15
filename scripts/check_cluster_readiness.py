#!/usr/bin/env python3
"""Report installation gaps without treating downloads as scientific qualification."""
import argparse
from datetime import datetime, timezone
import json
from pathlib import Path
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
    "mondo_disease_scope_registry",
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


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", type=Path, required=True)
    parser.add_argument("--state-dir", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    config = yaml.safe_load(args.config.read_text())
    resources = [{"key": key, "path": config.get(key), "status": inspect(config.get(key))} for key in RUNTIME_KEYS]
    if config.get("pext_enabled"):
        resources.append({"key": "pext_mean_bw", "path": config.get("pext_mean_bw"), "status": inspect(config.get("pext_mean_bw"))})
    vep_info = Path(config["vep_cache_dir"]) / "homo_sapiens_merged" / (str(config["vep_installed_version"]) + "_GRCh37") / "info.txt"
    resources.append({"key": "vep_cache_info", "path": str(vep_info), "status": inspect(vep_info)})
    receipts = {}
    for path in sorted(args.state_dir.glob("*.ready.json")):
        receipts[path.name] = json.loads(path.read_text())
    report = {
        "schema_version": 1, "checked_utc": datetime.now(timezone.utc).isoformat(),
        "status": "incomplete_installation", "config": str(args.config),
        "software_qualified": "environment.ready.json" in receipts,
        "resources": resources, "completed_stage_receipts": receipts,
        "blockers": [
            {"resource": "SpliceAI v1 + PrimateAI v0.2", "action": "Import the original licensed prescore files from user-provided paths or the intended BaseSpace account; do not drop their annotation/criteria."},
            {"resource": "MaveDB hg19", "action": "Discuss the unavailable assembly-specific annotation and each dependent calculation before omission."},
            {"resource": "Shared PriVA/CADD official VEP113", "action": "Qualify one annotation pass, transcript/peptide/GTF views and complete feature mapping. Packaged CADD VEP110 is reference-only."},
            {"resource": "Derived calibration caches", "action": "Complete gnomAD hg19 liftover, AlphaMissense/ClinVar statistics and other configured cache builders after inputs finish."},
            {"resource": "Family acceptance", "action": "User will provide family paths. Compare exact variant ordering and ACMG judgments with the same completed cache snapshot."},
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
