#!/usr/bin/env python3
"""Non-scientific compatibility glue: validate, isolate writes, invoke frozen PriVA.

No variant scoring, filtering, criteria, or sorting is implemented here. Existing
shell entry points remain authoritative until separately validated Rust ports.
"""
from __future__ import annotations

import argparse
import gzip
import hashlib
import json
import os
from pathlib import Path
import re
import shutil
import subprocess
import sys

import yaml

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
from spliceai_cache_input import inspect_header

ASSEMBLY = {"hg19": "hg19", "GRCh37": "hg19", "hg38": "hg38", "GRCh38": "hg38"}
PRIORITY_KEYS = (
    "clinvar_aa_stat", "clinvar_splice_stat", "all_intolerant_domains",
    "alphamissense_intolerant_domains", "pm1_regions_pkl", "alphamissense_vcf",
    "alphamissense_tranx_domain_map", "clinvar_patho_af_stat", "clinvar_patho_exon_af_stat",
    "interpro_mapping_pickle", "clingen_map", "gene_dosage_sensitivity",
    "hpo_condition_mechanism_json", "gene_nonlof_mechanism_json",
)
ANNOTATION_KEYS = (
    "ref_genome", "gnomad_vcf_chrX", "clinvar_vcf", "vep_cache_dir", "vep_plugins_dir",
    "vep_plugins_cachedir", "utr_annotator_file", "alphamissense_prescore",
    "spliceai_snv_prescore", "spliceai_indel_prescore", "primateai_prescore",
    "conservation_file", "splicevault_prescore", "loftee_repo", "human_ancestor_fasta",
    "loftee_conservation_file", "cadd_base_dir",
)
OPTIONAL_INPUTS = ("control_vcf", "alt_disease_vcf", "pp1_vcf", "pp1_ped",
                   "relevant_gene_list", "dispensable_gene_list", "cds_fasta_file")
SAFE_PATH = re.compile(r"^[A-Za-z0-9_./:+-]+$")
SAFE_ID = re.compile(r"^[A-Za-z0-9_][A-Za-z0-9_.-]*$")


class ContractError(ValueError):
    pass


class UniqueLoader(yaml.SafeLoader):
    pass


def unique_mapping(loader, node, deep=False):
    result = {}
    for key_node, value_node in node.value:
        key = loader.construct_object(key_node, deep=deep)
        if key in result:
            raise ContractError(f"Duplicate YAML key: {key}")
        result[key] = loader.construct_object(value_node, deep=deep)
    return result


UniqueLoader.add_constructor(yaml.resolver.BaseResolver.DEFAULT_MAPPING_TAG, unique_mapping)


def sha256(path):
    digest = hashlib.sha256()
    with open(path, "rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def json_digest(value):
    return hashlib.sha256(json.dumps(value, sort_keys=True, separators=(",", ":")).encode()).hexdigest()


def source_digest(root):
    root = Path(root)
    files = []
    for folder in ("scripts", "bin", "modules", "conf"):
        files.extend(path for path in (root / folder).rglob("*")
                     if path.is_file() and not path.name.endswith(".pyc") and "__pycache__" not in str(path))
    files.extend(root / name for name in ("main.nf", "nextflow.config"))
    digest = hashlib.sha256()
    # Match Nextflow's full path-string ordering. Path objects sort by parts,
    # which incorrectly puts conf/cadd/profile.json before conf/cadd.config.
    for path in sorted(files, key=lambda path: path.as_posix()):
        digest.update(str(path.relative_to(root)).encode())
        digest.update(b"\0")
        digest.update(path.read_bytes())
    return digest.hexdigest()


def software_readiness(mode):
    required = ["bash", "python", "python3", "bcftools", "tabix", "yq", "parallel", "tsv-pretty"]
    if mode == "raw":
        required += ["vep", "samtools", "conda"]
    missing = [tool for tool in required if shutil.which(tool) is None]
    if missing:
        raise ContractError("Runtime tools missing from PATH: " + ", ".join(missing))
    versions = {}
    for tool in required:
        entry = {"executable": str(Path(shutil.which(tool)).resolve())}
        if tool in {"python", "python3", "bcftools", "tabix", "yq", "samtools", "conda"}:
            result = subprocess.run([tool, "--version"], text=True, capture_output=True, timeout=30)
            entry["version_output"] = (result.stdout + result.stderr)[:4096]
            entry["version_exit_code"] = result.returncode
        versions[tool] = entry
    # The shell invokes python (not necessarily the adapter's python3), so check
    # the interpreter that will actually load the unchanged scientific modules.
    subprocess.run(["python", "-c", "import pandas,numpy,pysam,yaml,scipy"], check=True)
    return versions


def require_qualified_execution(mode, stub):
    if mode == "raw" and not stub:
        raise ContractError("Real raw annotation is gated: the unchanged CADD launcher may create Conda environments inside shared caches. Qualify enforced read-only cache/container execution before enabling this entry point. The first acceptance run can use annotated mode with the original paired CADD TSV.")


def read_config(path):
    with open(path) as handle:
        config = yaml.load(handle, Loader=UniqueLoader)
    if not isinstance(config, dict):
        raise ContractError("Legacy config must be a YAML mapping")
    if config.get("assembly") not in ASSEMBLY:
        raise ContractError("Config assembly must be hg19/hg38/GRCh37/GRCh38")
    # Annotated imports skip the annotation shell. Inspect the actual resources
    # at the original two keys here too, including cached launches and tasks;
    # an installer setting cannot confer scientific qualification.
    for key in ("spliceai_snv_prescore", "spliceai_indel_prescore"):
        value = config.get(key)
        if not value or not isinstance(value, str) or not Path(value).is_file():
            continue  # Keep existing missing-resource diagnostics at their stage.
        try:
            detected = inspect_header(value)
        except (OSError, EOFError, ValueError) as error:
            raise ContractError(f"Cannot verify {key} header: {error}") from error
        if detected["format"] != "vcf_v1":
            raise ContractError("SpliceAI2 is prepared for annotation inspection only; legacy PriVA ACMG/ClinVar criteria are not qualified for v2. Use original v1 resources for this scientific reference.")
    return config


def checked_path(value, key, absolute=True):
    if not isinstance(value, str) or not value or not SAFE_PATH.fullmatch(value):
        raise ContractError(f"{key}: required path is empty or contains whitespace/shell metacharacters")
    path = Path(value)
    if absolute and not path.is_absolute():
        raise ContractError(f"{key}: expected an absolute path, got {value}")
    if not path.exists():
        raise ContractError(f"{key}: missing path {value}")
    if path.is_file():
        with path.open("rb") as handle:
            if handle.read(200).startswith(b"version https://git-lfs.github.com/spec/v1"):
                raise ContractError(f"{key}: Git LFS pointer must be materialized: {value}")
    return path


def runtime_resources(config, root, mode, vcf_hub_mode="off", cadd_hub_mode="copy"):
    """Mirror actual shell consumers, including their repository-relative files."""
    paths = {}
    for key in PRIORITY_KEYS + (ANNOTATION_KEYS if mode == "raw" else ()):
        paths[key] = checked_path(config.get(key), key)
    for key in OPTIONAL_INPUTS:
        if config.get(key):
            paths[key] = checked_path(config[key], key)
    if bool(config.get("pp1_vcf")) != bool(config.get("pp1_ped")):
        raise ContractError("pp1_vcf and pp1_ped must be supplied together")
    fixed = {
        "hpo_assertions_used_by_shell": root / "data/hpo/genes_to_phenotype.assertions.tsv.gz",
        "alphamissense_mean_used_by_shell": root / "data/alphamissense/alphamissense_mean_score.tsv",
        "clinvar_gene_stat_used_by_shell": root / "data/ClinVar" / Path(str(config.get("clinvar_gene_stat", ""))).name,
        "repeat_regions_used_by_shell": root / "data/repeats" / Path(str(config.get("repeat_region_file_name", ""))).name,
        "hgnc_used_by_scientific_code": root / "data/hgnc/non_alt_loci_set.tsv",
    }
    if not config.get("dispensable_gene_list"):
        fixed["dispensable_genes_default"] = root / "data/dispensable_genes/dispensable_gene_list.txt"
    if not config.get("clinvar_gene_stat") or not config.get("repeat_region_file_name"):
        raise ContractError("clinvar_gene_stat and repeat_region_file_name are required")
    if mode == "raw":
        fixed["primary_contigs"] = root / "data/liftover/ucsc_GRC.primary.contigs.tsv"
        fixed["contig_name_mapping"] = root / "data/liftover/ucsc_to_GRC.contig.map.tsv"
        if not config.get("loeuf_prescore"):
            raise ContractError("loeuf_prescore is required")
        fixed["loeuf_used_by_shell"] = root / "data/loeuf" / Path(config["loeuf_prescore"]).name
        # gnomAD is read chromosome-by-chromosome from the chrX template.
        template = str(paths["gnomad_vcf_chrX"])
        if "chrX" not in template:
            raise ContractError("gnomad_vcf_chrX must contain chrX for the legacy template")
        for chromosome in list(map(str, range(1, 23))) + ["Y"]:
            key = f"gnomad_chr{chromosome}"
            paths[key] = checked_path(template.replace("chrX", f"chr{chromosome}"), key)
        fixed["reference_fai"] = Path(str(paths["ref_genome"]) + ".fai")
        if ASSEMBLY[config["assembly"]] == "hg38":
            for key in ("gerp_bigwig", "mavedb_file"):
                paths[key] = checked_path(config.get(key), key)
        if config.get("pext_enabled") is True:
            paths["pext_mean_bw"] = checked_path(config.get("pext_mean_bw"), "pext_mean_bw")
            if config.get("pext_tissues"):
                paths["pext_tissue_bw_dir"] = checked_path(config.get("pext_tissue_bw_dir"), "pext_tissue_bw_dir")
        for key, mode_value in (("hub_vcf_file", vcf_hub_mode), ("hub_cadd_file", cadd_hub_mode)):
            if mode_value != "off" and config.get(key):
                path = Path(config[key])
                if not path.is_absolute():
                    path = Path(config.get("hub_file_dir", "")) / path
                paths[key] = checked_path(str(path), key)
    if ASSEMBLY[config["assembly"]] == "hg38":
        fixed["mavedb_assay_metadata_used_by_shell"] = root / "data/MaveDB/assay_mapping_result.neutral_bands.tsv"
    for key, path in fixed.items():
        paths[key] = checked_path(str(path), key)
    if "cds_fasta_file" in paths:
        # Copy this mutable lazy cache into each prioritization task; don't write
        # a newly parsed pickle beside the shared FASTA.
        cache = Path(str(paths["cds_fasta_file"]).replace(".fa.gz", ".cds_cache.pkl.gz"))
        if cache != paths["cds_fasta_file"] and cache.exists():
            paths["cds_parsed_cache"] = checked_path(str(cache), "cds_parsed_cache")
    for key, path in list(paths.items()):
        if str(path).endswith((".vcf.gz", ".vcf.bgz")):
            # Legacy tools rebuild missing tabix indexes; require a ready index.
            paths[key + "_index"] = checked_path(str(path) + ".tbi", key + "_index")
    return paths


def inventory_record(path, checksum=False):
    path = Path(path)
    stat = path.stat()
    record = {"path": str(path), "resolved_path": str(path.resolve()),
              "type": "directory" if path.is_dir() else "file", "size": stat.st_size,
              "mtime_ns": stat.st_mtime_ns}
    if checksum and path.is_file():
        record["sha256"] = sha256(path)
    return record


def inventory(paths, checksum=False):
    records = {}
    for key, path in sorted(paths.items()):
        records[key] = inventory_record(path, checksum)
        if path.is_dir():
            # Include every dependency below directory resources without hashing
            # TB of bytes by default. Symlinks are recorded and their target stat
            # included, but directory symlinks are not silently traversed.
            for child in sorted(path.rglob("*")):
                if child.is_symlink() and child.is_dir():
                    raise ContractError(f"Directory symlink in resource tree needs explicit resolved config: {child}")
                if child.is_file():
                    records[f"{key}/{child.relative_to(path)}"] = inventory_record(child, checksum)
    return records


def write_json(path, value):
    Path(path).write_text(json.dumps(value, indent=2, sort_keys=True) + "\n")


def build_manifest(args):
    config = read_config(args.config)
    root = checked_path(args.legacy_root, "legacy_root").resolve()
    paths = runtime_resources(config, root, args.mode, args.vcf_hub_mode, args.cadd_hub_mode)
    payload = {"schema_version": 1, "assembly": ASSEMBLY[config["assembly"]], "mode": args.mode,
               "config_sha256": sha256(args.config), "legacy_root": str(root),
               "identity_kind": "sha256" if args.checksum else "size-mtime-inventory",
               "vcf_hub_mode": args.vcf_hub_mode, "cadd_hub_mode": args.cadd_hub_mode,
               "resources": inventory(paths, args.checksum)}
    payload["bundle_id"] = json_digest(payload)
    write_json(args.output, payload)


def verify_manifest(manifest, config_path, root, mode, vcf_hub_mode, cadd_hub_mode):
    config = read_config(config_path)
    if manifest.get("schema_version") != 1 or not manifest.get("resources"):
        raise ContractError("Resource manifest is incomplete; create it with priva_task.py manifest")
    recorded_id = manifest.get("bundle_id")
    if recorded_id != json_digest({k: v for k, v in manifest.items() if k != "bundle_id"}):
        raise ContractError("Resource manifest bundle_id does not match its content")
    expected = {"config_sha256": sha256(config_path), "legacy_root": str(root), "mode": mode,
                "assembly": ASSEMBLY[config["assembly"]], "vcf_hub_mode": vcf_hub_mode,
                "cadd_hub_mode": cadd_hub_mode}
    for key, value in expected.items():
        if manifest.get(key) != value:
            raise ContractError(f"Resource manifest mismatch: {key}")
    paths = runtime_resources(config, root, mode, vcf_hub_mode, cadd_hub_mode)
    current = inventory(paths, checksum=False)
    if set(current) != set(manifest["resources"]):
        raise ContractError("Resource inventory membership changed; freeze and rebuild manifest")
    for key, record in current.items():
        expected_record = {k: v for k, v in manifest["resources"][key].items() if k != "sha256"}
        if record != expected_record:
            raise ContractError(f"Resource changed since manifest creation: {key}")
    return paths


def vcf_samples(path):
    with gzip.open(path, "rt") as handle:
        for line in handle:
            if line.startswith("##SpliceAI2") or line.rstrip("\r\n") == "##source=PriVA_SpliceAI2_projection":
                raise ContractError("SpliceAI2-marked input is not qualified for legacy PriVA criteria")
            if line.startswith("#CHROM\t"):
                fields = line.rstrip("\n").split("\t")
                if len(fields) < 10 or fields[8] != "FORMAT":
                    raise ContractError("Family input VCF must have FORMAT and sample columns")
                return fields[9:]
            if not line.startswith("#"):
                break
    raise ContractError("Input VCF has no #CHROM header")


def validate_ped(path, family, samples):
    rows = [line.split("\t") for line in Path(path).read_text().splitlines()
            if line.strip() and not line.lstrip().startswith("#")]
    if not rows or any(len(row) < 6 for row in rows):
        raise ContractError("Legacy Python requires a tab-delimited PED with at least six columns and one person")
    selected = [row for row in rows if row[0] == family]
    if not selected:
        raise ContractError("The selected family is absent from the PED")
    ids = [row[1] for row in rows]
    if len(ids) != len(set(ids)) or len(samples) != len(set(samples)):
        raise ContractError("Duplicate sample identifier in PED or VCF")
    if any(not SAFE_ID.fullmatch(value) for value in ids):
        raise ContractError("PED sample identifier contains unsupported shell characters")
    if set(ids) - set(samples):
        raise ContractError("Every PED sample must be present in the VCF")
    if any(row[4] not in {"0", "1", "2"} or row[5] not in {"1", "2"} for row in rows):
        raise ContractError("PED requires sex 0/1/2 and known affected status 1/2")
    if selected[0][5] != "2":
        raise ContractError("Legacy PriVA requires the affected proband first; PED order is never rewritten")
    return [row[1] for row in selected]


def validate_cadd(path):
    with open(path) as handle:
        fields = handle.readline().rstrip("\n").split("\t")
    required = {"#Chrom", "Pos", "Ref", "Alt", "PHRED", "FeatureID"}
    if not required.issubset(fields) or len(fields) != len(set(fields)):
        raise ContractError("CADD first row must have unique #Chrom, Pos, Ref, Alt, PHRED, FeatureID columns; remove no rows implicitly")


def validate(args):
    require_qualified_execution(args.mode, args.stub)
    config = read_config(args.config)
    for value in (args.family, args.base):
        if not SAFE_ID.fullmatch(value):
            raise ContractError("Unsafe family/base identifier")
    if ASSEMBLY[config["assembly"]] != ASSEMBLY[args.assembly]:
        raise ContractError("Requested assembly disagrees with legacy config")
    if args.mode == "annotated" and not args.cadd:
        raise ContractError("Annotated mode requires the paired CADD TSV")
    samples = vcf_samples(args.vcf)
    family_samples = validate_ped(args.ped, args.family, samples)
    if not Path(args.vcf_index).is_file():
        raise ContractError("Input VCF index is missing")
    if args.cadd:
        validate_cadd(args.cadd)
    manifest = json.loads(Path(args.resources).read_text())
    root = Path(args.legacy_root).resolve()
    if source_digest(root) != args.source_digest:
        raise ContractError("Source changed after Nextflow launch; freeze the checkout before retrying")
    versions = {}
    if not args.stub:
        verify_manifest(manifest, args.config, root, args.mode, args.vcf_hub_mode, args.cadd_hub_mode)
        versions = software_readiness(args.mode)
        subprocess.run(["bcftools", "index", "--nrecords", args.vcf], check=True, stdout=subprocess.DEVNULL)
    bundle = Path("bundle")
    bundle.mkdir(exist_ok=False)
    for source, name in ((args.config, "legacy.yaml"), (args.ped, "family.ped"), (args.resources, "resources.json")):
        shutil.copy2(source, bundle / name)
    contract = {"schema_version": 1, "family": args.family, "base": args.base,
                "assembly": ASSEMBLY[args.assembly], "mode": args.mode, "legacy_root": str(root),
                "source_digest": args.source_digest, "software_id": args.software_id,
                "resource_bundle_id": manifest.get("bundle_id", "stub"),
                "vcf_hub_mode": args.vcf_hub_mode, "cadd_hub_mode": args.cadd_hub_mode,
                "samples": samples, "family_samples": family_samples, "stub": args.stub,
                "input_sha256": sha256(args.vcf), "ped_sha256": sha256(args.ped),
                "provided_cadd_sha256": sha256(args.cadd) if args.cadd else None}
    write_json(bundle / "contract.json", contract)
    write_json(bundle / "software-observed.json", versions)


def launch(args):
    """Runs in the driver on every launch, including a completely cached resume."""
    require_qualified_execution(args.mode, args.stub)
    root = Path(args.legacy_root).resolve()
    if source_digest(root) != args.source_digest:
        raise ContractError("Source changed while launching")
    config = read_config(args.config)
    if ASSEMBLY[config["assembly"]] != ASSEMBLY[args.assembly]:
        raise ContractError("Requested assembly disagrees with legacy config")
    manifest = json.loads(Path(args.resources).read_text())
    if not args.stub:
        verify_manifest(manifest, args.config, root, args.mode, args.vcf_hub_mode, args.cadd_hub_mode)
    fingerprint = {"schema_version": 1, "source_digest": args.source_digest, "software_id": args.software_id,
                   "family": args.family, "assembly": ASSEMBLY[args.assembly], "mode": args.mode,
                   "config_sha256": sha256(args.config), "ped_sha256": sha256(args.ped),
                   "manifest_sha256": sha256(args.resources), "input": inventory_record(Path(args.vcf).resolve()),
                   "input_index": inventory_record(Path(args.vcf_index).resolve()),
                   "provided_cadd": inventory_record(Path(args.cadd).resolve()) if args.cadd else None,
                   "vcf_hub_mode": args.vcf_hub_mode, "cadd_hub_mode": args.cadd_hub_mode,
                   "annotation_cpus": args.annotation_cpus, "family_cpus": args.family_cpus, "stub": args.stub}
    marker = {"run_identity": json_digest(fingerprint), "contract": fingerprint}
    outdir = Path(args.outdir).resolve()
    outdir.mkdir(parents=True, exist_ok=True)
    path = outdir / ".priva-run.json"
    if path.exists():
        if json.loads(path.read_text()) != marker:
            raise ContractError("Output directory belongs to different inputs/source/resources/settings; use a new --outdir")
    else:
        # Nextflow may create trace/report directories before entering workflow.
        if any(child.name != "pipeline_info" for child in outdir.iterdir()):
            raise ContractError("Output directory is not empty and has no matching PriVA run identity")
        with path.open("x") as handle:
            json.dump(marker, handle, indent=2, sort_keys=True)
            handle.write("\n")


def copy_writable(source, destination):
    destination = Path(destination)
    destination.parent.mkdir(parents=True, exist_ok=True)
    shutil.copy2(source, destination)
    destination.chmod(destination.stat().st_mode | 0o200)
    return destination.absolute()


def private_config(bundle, contract, phase, cpus):
    config = read_config(bundle / "legacy.yaml")
    result = Path("result").absolute()
    result.mkdir(exist_ok=True)
    temp = Path("tmp").absolute()
    temp.mkdir(exist_ok=True)
    config.update(input_vcf=str(result / (contract["base"] + ".vcf.gz")),
                  output_dir=str(result), tmp_dir=str(temp),
                  ped_file=str(copy_writable(bundle / "family.ped", result / "family.ped")),
                  threads=cpus, threads_per_fam=cpus, assembly=contract["assembly"],
                  hub_vcf_file="", hub_cadd_file="", hub_file_dir=str(result / "private_hubs"))
    for key in OPTIONAL_INPUTS:
        source = config.get(key)
        if not source:
            continue
        if phase != "prioritization" and key == "cds_fasta_file":
            continue
        dest = copy_writable(source, result / "private_inputs" / key / Path(source).name)
        config[key] = str(dest)
        if str(source).endswith((".vcf.gz", ".vcf.bgz")):
            copy_writable(str(source) + ".tbi", str(dest) + ".tbi")
        if key == "cds_fasta_file":
            original_cache = Path(source.replace(".fa.gz", ".cds_cache.pkl.gz"))
            if original_cache != Path(source) and original_cache.exists():
                copy_writable(original_cache, str(dest).replace(".fa.gz", ".cds_cache.pkl.gz"))
    if phase == "annotation" and contract["mode"] == "raw":
        original = read_config(bundle / "legacy.yaml")
        for key, option in (("hub_vcf_file", "vcf_hub_mode"), ("hub_cadd_file", "cadd_hub_mode")):
            if contract[option] == "off" or not original.get(key):
                continue
            source = Path(original[key])
            if not source.is_absolute():
                source = Path(original["hub_file_dir"]) / source
            destination = copy_writable(source, result / "private_hubs" / source.name)
            config[key] = str(destination)
            if key == "hub_vcf_file":
                copy_writable(str(source) + ".tbi", str(destination) + ".tbi")
    path = Path("task.yaml").absolute()
    path.write_text(yaml.safe_dump(config, sort_keys=False))
    return path, temp


def require_outputs(paths):
    for path in paths:
        if not Path(path).is_file():
            raise ContractError(f"Legacy task returned without required output: {path}")


def run(args):
    bundle = Path(args.bundle)
    contract = json.loads((bundle / "contract.json").read_text())
    require_qualified_execution(contract["mode"], args.stub)
    if bool(contract["stub"]) != bool(args.stub):
        raise ContractError("Stub and real task contracts cannot be mixed")
    if args.cpus < 1:
        raise ContractError("cpus must be positive")
    root = Path(contract["legacy_root"])
    if source_digest(root) != contract["source_digest"]:
        raise ContractError("Scientific or adapter source changed after validation; task refused")
    if not args.stub:
        manifest = json.loads((bundle / "resources.json").read_text())
        verify_manifest(manifest, bundle / "legacy.yaml", root, contract["mode"],
                        contract["vcf_hub_mode"], contract["cadd_hub_mode"])
    result = Path("result").absolute()
    result.mkdir(exist_ok=False)
    base = contract["base"]
    family_base = f"{base}.anno.{contract['family']}.filtered"
    names = {"annotation": base + (".anno" if contract["mode"] == "annotated" else ""),
             "filtration": base + ".anno", "prioritization": family_base}
    input_vcf = copy_writable(args.vcf, result / (names[args.phase] + ".vcf.gz"))
    copy_writable(args.vcf_index, str(input_vcf) + ".tbi")
    if args.cadd:
        cadd_base = family_base if args.phase == "prioritization" else base + ".anno"
        cadd_path = copy_writable(args.cadd, result / (cadd_base + ".cadd.tsv"))
    elif args.phase != "annotation":
        raise ContractError("Filtration/prioritization require the paired CADD TSV")
    outputs = []
    if args.stub:
        # Tiny syntactic fixtures; no scientific acceptance is claimed.
        if args.phase == "annotation":
            destination = result / (base + ".anno.vcf.gz")
            if input_vcf != destination:
                copy_writable(input_vcf, destination)
                copy_writable(str(input_vcf) + ".tbi", str(destination) + ".tbi")
            (result / (base + ".anno.cadd.tsv")).write_text("#Chrom\tPos\tRef\tAlt\tRawScore\tPHRED\tFeatureID\n1\t10\tA\tG\t1.0\t10\tSTUB\n")
        elif args.phase == "filtration":
            destination = result / (family_base + ".vcf.gz")
            copy_writable(input_vcf, destination)
            copy_writable(str(input_vcf) + ".tbi", str(destination) + ".tbi")
            copy_writable(cadd_path, result / (family_base + ".cadd.tsv"))
        else:
            (result / (family_base + ".tsv")).write_text("Chr\tStart\tRef\tAlt\tFeature\tsort_index\tvariant_rank\tACMG_class\tACMG_criteria\n1\t10\tA\tG\tSTUB\t1\t1\tSTUB_ONLY\tSTUB_ONLY\n")
            (result / (family_base + ".acmg.tsv")).write_text("variant_id\tSTUB_ONLY\n1:10:A:G\t1\n")
    else:
        config_path, temp = private_config(bundle, contract, args.phase, args.cpus)
        env = os.environ.copy()
        env.update(TMPDIR=str(temp), TMP=str(temp), TEMP=str(temp), PYTHONDONTWRITEBYTECODE="1",
                   OMP_NUM_THREADS=str(args.cpus), OPENBLAS_NUM_THREADS=str(args.cpus),
                   MKL_NUM_THREADS=str(args.cpus), NUMEXPR_NUM_THREADS=str(args.cpus))
        if args.phase == "annotation" and contract["mode"] == "raw":
            command = ["bash", str(root / "scripts/annotation_vcf.sh"), "main_workflow", "-i", str(input_vcf), "--config", str(config_path)]
        elif args.phase == "filtration":
            command = ["bash", str(root / "scripts/filtration_vcf_per_fam.sh"), "-v", str(input_vcf), "-f", contract["family"], "-c", str(config_path)]
        elif args.phase == "prioritization":
            command = ["bash", str(root / "scripts/prioritization_vcf_per_fam.sh"), "-i", str(input_vcf), "-f", contract["family"], "-c", str(config_path), "-t", str(cadd_path)]
        else:
            command = None  # Explicit preannotated import retains input bytes.
        if command:
            subprocess.run(command, check=True, env=env)
    if args.phase == "annotation":
        outputs = [result / (base + suffix) for suffix in (".anno.vcf.gz", ".anno.vcf.gz.tbi", ".anno.cadd.tsv")]
    elif args.phase == "filtration":
        outputs = [result / (family_base + suffix) for suffix in (".vcf.gz", ".vcf.gz.tbi", ".cadd.tsv")]
    else:
        outputs = [result / (family_base + suffix) for suffix in (".tsv", ".acmg.tsv")]
    require_outputs(outputs)
    if not args.stub:
        verify_manifest(manifest, bundle / "legacy.yaml", root, contract["mode"],
                        contract["vcf_hub_mode"], contract["cadd_hub_mode"])
    write_json(f"{args.phase}.provenance.json", {**contract, "phase": args.phase, "cpus": args.cpus,
        "outputs": [{"name": path.name, "sha256": sha256(path), "size": path.stat().st_size} for path in outputs],
        "cache_isolation": "shared resources read; mutable VCF/CADD hubs, optional cohort VCFs and CDS lazy cache copied into task"})


def parser():
    command = argparse.ArgumentParser(description=__doc__)
    sub = command.add_subparsers(dest="command", required=True)
    manifest = sub.add_parser("manifest", help="Inventory existing prepared caches; never installs or rewrites them")
    manifest.add_argument("--config", required=True)
    manifest.add_argument("--legacy-root", required=True)
    manifest.add_argument("--mode", choices=["raw", "annotated"], default="raw")
    manifest.add_argument("--vcf-hub-mode", choices=["off", "qualified-copy"], default="off")
    manifest.add_argument("--cadd-hub-mode", choices=["off", "copy"], default="copy")
    manifest.add_argument("--checksum", action="store_true", help="Read/hash every file; potentially TB of I/O")
    manifest.add_argument("--output", required=True)
    manifest.set_defaults(func=build_manifest)
    validate_cmd = sub.add_parser("validate")
    for key in ("vcf", "vcf-index", "ped", "config", "resources", "family", "base", "assembly", "mode", "legacy-root", "source-digest", "software-id", "vcf-hub-mode", "cadd-hub-mode"):
        validate_cmd.add_argument("--" + key, required=True)
    validate_cmd.add_argument("--cadd")
    validate_cmd.add_argument("--stub", action="store_true")
    validate_cmd.set_defaults(func=validate)
    launch_cmd = sub.add_parser("launch", help="Uncached launch-time resource verification and output identity guard")
    for key in ("vcf", "vcf-index", "ped", "config", "resources", "family", "assembly", "mode", "legacy-root", "source-digest", "software-id", "vcf-hub-mode", "cadd-hub-mode", "outdir"):
        launch_cmd.add_argument("--" + key, required=True)
    launch_cmd.add_argument("--annotation-cpus", type=int, required=True)
    launch_cmd.add_argument("--family-cpus", type=int, required=True)
    launch_cmd.add_argument("--cadd")
    launch_cmd.add_argument("--stub", action="store_true")
    launch_cmd.set_defaults(func=launch)
    run_cmd = sub.add_parser("run")
    run_cmd.add_argument("--phase", choices=["annotation", "filtration", "prioritization"], required=True)
    for key in ("bundle", "vcf", "vcf-index"):
        run_cmd.add_argument("--" + key, required=True)
    run_cmd.add_argument("--cadd")
    run_cmd.add_argument("--cpus", type=int, required=True)
    run_cmd.add_argument("--stub", action="store_true")
    run_cmd.set_defaults(func=run)
    return command


if __name__ == "__main__":
    try:
        arguments = parser().parse_args()
        arguments.func(arguments)
    except (ContractError, OSError, ValueError, subprocess.CalledProcessError) as exc:
        print(f"PriVA compatibility contract failed: {exc}", file=sys.stderr)
        sys.exit(2)
