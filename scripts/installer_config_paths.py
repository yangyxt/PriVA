#!/usr/bin/env python3
"""Plan managed cache destinations and adopt explicit source bytes before selecting them."""
from __future__ import annotations

import argparse
from datetime import datetime, timezone
import fcntl
import hashlib
import importlib.util
import json
import os
from pathlib import Path
import re
import shutil
import subprocess
import sys
import tempfile

import yaml


FOLDERS = {
    "vep_cache": ("VEP_cache", "vep_cache_dir"),
    "vep_plugin": ("VEP_plugin", "vep_plugins_dir"),
    "vep_plugin_cache": ("VEP_plugin_cache", "vep_plugins_cachedir"),
    "gnomad": ("gnomAD", "gnomad_vcf_dir"),
    "clinvar": ("ClinVar", "clinvar_vcf_dir"),
    "cadd": ("CADD", "cadd_parent_dir"),
}
BUNDLED_PATHS = {
    "conda_env_yaml": "acmg_conda_rough.yml",
    "hgnc_table": "data/hgnc/hgnc_non_alt_loci_set_sha256_890b18743cc8.tsv",
    "dispensable_gene_list": "data/dispensable_genes/dispensable_gene_list.txt",
    "gene_mechanism_builder_script": "scripts/build_gene_pathogenic_mechanism_cache.py",
    "gene_nonlof_mechanism_builder_script": "scripts/build_gene_nonlof_mechanism_cache.py",
    "gene_nonlof_mechanism_schema": "data/patho_mechanism/gene_nonlof_mechanism_curated_assertions.schema.json",
    "hpo_condition_mechanism_builder_script": "scripts/build_hpo_condition_mechanism_cache.py",
    "gofcards_clinvar_injector_script": "scripts/inject_clinvar_into_gofcards.py",
    "gofcards_builder_script": "scripts/build_gofcards_exact_gof_cache.py",
    "gofcards_mechanism_review_tsv": "data/gofcards/gofcards_mechanism_reviews.tsv",
    "hpo_assertion_builder_script": "scripts/collapse_HPO_anno.py",
    "mondo_scope_builder_script": "scripts/build_mondo_disease_scope.py",
    "mondo_disease_scope_overrides": "data/mondo/disease_scope_overrides.tsv",
}
EMPTY = (None, "", "null")
BUNDLED_CACHES = {
    "gofcards_exact_gof_cache": "data/gofcards/gofcards_exact_gof.json.gz",
    "gene_nonlof_mechanism_json": "data/patho_mechanism/gene_nonlof_mechanism_curated_assertions.json.gz",
    "hpo_condition_mechanism_json": "data/patho_mechanism/hpo_condition_mechanism_cache.json.gz",
    "ddg2p_mechanism_evidence": "data/patho_mechanism/gene_pathogenic_mechanism_evidence.tsv",
}
REPLACEMENT_TABLE_KEYS = {
    "hgnc_table", "dispensable_gene_list", "gofcards_mechanism_review_tsv",
    "mondo_disease_scope_overrides", "gene_nonlof_mechanism_schema",
}
FORCE_KEYS = {
    "gofcards_exact_gof_cache": ("PRIVA_FORCE_GOFCARDS_CACHE", "PRIVA_FORCE_GOFCARDS_CLINVAR_INJECTION"),
    "gene_nonlof_mechanism_json": ("PRIVA_FORCE_NONLOF_MECHANISM_CACHE",),
    "hpo_condition_mechanism_json": ("PRIVA_FORCE_HPO_MECHANISM_JSON",),
    "ddg2p_mechanism_evidence": ("PRIVA_FORCE_GENE_MECHANISM_CACHE",),
}

# These are installation resources, not workflow inputs/output directories.
# Existing derived files retain their bytes and scientific companion receipts;
# missing derived destinations are built by their existing scientific procedures.
GENERATED_CACHE_KEYS = {
    "clinvar_vcf", "clinvar_aa_stat", "clinvar_splice_stat", "clinvar_patho_af_stat",
    "clinvar_patho_exon_af_stat", "clinvar_gene_stat", "alphamissense_vcf",
    "alphamissense_vep_vcf", "alphamissense_vep_base_vcf", "alphamissense_pd_stat", "alphamissense_tranx_domain_map",
    "alphamissense_intolerant_domains", "alphamissense_intolerant_motifs",
    "all_intolerant_domains", "pm1_regions_pkl", "interpro_mapping_pickle",
    "clingen_map", "gene_dosage_sensitivity", "gene_nonlof_mechanism_json",
    "hpo_condition_mechanism_json", "ddg2p_mechanism_evidence", "gofcards_exact_gof_cache",
    "hpo_assertions", "mondo_disease_scope_registry", "clinvar_vcv_match_audit",
}
DIRECTORY_CACHE_KEYS = {
    *(key for _, key in FOLDERS.values()), "loftee_parent_dir", "loftee_repo", "cadd_base_dir",
    "gene_mechanism_cache_dir", "gene_mechanism_raw_dir", "gofcards_workdir", "hpo_raw_dir",
    "pext_tissue_bw_dir",
}
GENERATED_DIRECTORY_KEYS = {
    *(key for _, key in FOLDERS.values()), "loftee_parent_dir", "gene_mechanism_cache_dir",
    "gene_mechanism_raw_dir", "gofcards_workdir", "hpo_raw_dir",
}
AUTO_SOURCE_KEYS = {
    "ref_genome", "gnomad_vcf_chrX", "alphamissense_prescore", "utr_annotator_file",
    "human_ancestor_fasta", "loftee_conservation_file", "gerp_bigwig", "conservation_file",
    "splicevault_prescore", "cds_fasta_file", "mondo_obo", "clinvar_vcv_xml", "mavedb_file",
    "pext_mean_bw", "pext_tissue_bw_dir", "loftee_repo", "cadd_base_dir",
}
GOFCARDS_BUILD_ONLY_KEYS = {
    "gofcards_mechanism_review_tsv",
}
# Explicit, known scientific companions only. Contents are copied unchanged;
# copying is not a new assertion that the scientific cache is qualified.
FILE_COMPANIONS = (".tbi", ".csi", ".fai", ".gzi", ".provenance.json", ".nm_to_enst.json",
                   ".ready.json", ".source.json", ".annotation-base.json")


def validate_repeat_filename(value) -> None:
    """The existing runtime accepts a filename beneath data/repeats only."""
    if value not in EMPTY and (not isinstance(value, str) or value in (".", "..")
                              or not re.fullmatch(r"[^/\\\r\n]+", value)):
        raise ValueError("repeat_region_file_name must be one filename beneath data/repeats, not an absolute or relative path")


def reverse_chain_source(value, source_root: Path | None) -> Path:
    """Keep the legacy bundled filename or an explicit absolute source path."""
    root = source_root or Path(__file__).resolve().parents[1]
    if value in EMPTY:
        return root / "data/liftover/hg38ToHg19.over.chain.gz"
    if not isinstance(value, str):
        raise ValueError("hg38_hg19_chain must be a bundled filename or an absolute supplied path")
    selected = Path(value)
    if selected.is_absolute():
        return selected
    if not re.fullmatch(r"[A-Za-z0-9_][A-Za-z0-9_.-]*", value):
        raise ValueError("hg38_hg19_chain must be a bundled filename or an absolute supplied path")
    return root / "data/liftover" / value


def resource_groups(config: dict, directories: dict[str, str]) -> dict[str, Path]:
    root = checked_storage(config["cache_root"], "cache_root")
    plugin = Path(directories["vep_plugins_cachedir"])
    result = {}
    groups = {
        plugin / "AlphaMissense": ("alphamissense_prescore", "alphamissense_vcf", "alphamissense_vep_vcf", "alphamissense_vep_base_vcf",
            "alphamissense_pd_stat", "alphamissense_tranx_domain_map", "alphamissense_intolerant_domains",
            "alphamissense_intolerant_motifs", "all_intolerant_domains", "pm1_regions_pkl"),
        Path(directories["clinvar_vcf_dir"]): ("clinvar_vcf", "clinvar_aa_stat", "clinvar_splice_stat",
            "clinvar_patho_af_stat", "clinvar_patho_exon_af_stat", "clinvar_gene_stat"),
        Path(directories["gnomad_vcf_dir"]): ("gnomad_vcf_chrX",),
        plugin / "UTRannotator": ("utr_annotator_file",),
        plugin / "LoFtee": ("human_ancestor_fasta", "loftee_conservation_file", "gerp_bigwig"),
        Path(directories["vep_plugins_dir"]) / "LoFtee_repos": ("loftee_parent_dir", "loftee_repo"),
        plugin / "Conservation": ("conservation_file",),
        plugin / "SpliceVault": ("splicevault_prescore",),
        plugin / "PrimateAI": ("primateai_prescore",),
        plugin / "SpliceAI": ("spliceai_snv_prescore", "spliceai_indel_prescore"),
        plugin / "Ensembl_CDS": ("cds_fasta_file",),
        plugin / "MaveDB": ("mavedb_file",),
        plugin / "pext": ("pext_mean_bw", "pext_tissue_bw_dir"),
        Path(directories["cadd_parent_dir"]): ("cadd_base_dir",),
        root / "InterPro": ("interpro_mapping_pickle",),
        root / "ClinGen": ("clingen_map", "gene_dosage_sensitivity"),
        root / "LOEUF": ("loeuf_prescore",),
        root / "reference": ("ref_genome", "gofcards_hg19_fasta", "gofcards_hg38_fasta"),
        root / "liftOver": ("gofcards_hg19_to_hg38_chain", "hg38_hg19_chain"),
        root / "patho_mechanism": ("gene_mechanism_cache_dir", "gene_mechanism_raw_dir",
            "gene_nonlof_mechanism_json", "hpo_condition_mechanism_json", "ddg2p_mechanism_evidence",
            "clinvar_vcv_match_audit"),
        root / "gofcards": ("gofcards_workdir", "gofcards_exact_gof_cache"),
        root / "patho_mechanism/raw/clinvar_vcv": ("clinvar_vcv_xml",),
        root / "hpo": ("hpo_raw_dir", "hpo_assertions"),
        root / "mondo": ("mondo_obo", "mondo_disease_scope_registry"),
    }
    for folder, keys in groups.items():
        result.update({key: folder for key in keys})
    for key in REPLACEMENT_TABLE_KEYS - {"hgnc_table"}:
        result[key] = root / Path(BUNDLED_PATHS[key]).parent
    return result


def bundled_resource(key: str, value: str, source_root: Path | None) -> Path | None:
    """Recognize only canonical checkout snapshots, using existence/size alone."""
    if source_root is None:
        return None
    relative = dict(BUNDLED_CACHES)
    relative.update({key: BUNDLED_PATHS[key] for key in REPLACEMENT_TABLE_KEYS})
    relative.update({"interpro_mapping_pickle": "data/InterPro/Interpro_entry_mapping.pkl.gz",
                     "clinvar_gene_stat": "data/ClinVar/clinvar_2star_stats.pkl.gz",
                     "hg38_hg19_chain": "data/liftover/hg38ToHg19.over.chain.gz"})
    if key == "hg38_hg19_chain":
        value = str(reverse_chain_source(value, source_root))
    options = [relative[key]] if key in relative else []
    if key == "hgnc_table":
        options.append("data/hgnc/non_alt_loci_set.tsv")
    if key == "loeuf_prescore":
        options = ["data/loeuf/loeuf_dataset.tsv.gz", "data/loeuf/loeuf_dataset_hg38.tsv.gz"]
    if key in ("clingen_map", "gene_dosage_sensitivity"):
        stem, suffix = ("clingen_map", "pkl.gz") if key == "clingen_map" else ("gene_dosage_sensitivity", "tsv")
        options = [f"data/clingen/{stem}.{assembly}.{suffix}" for assembly in ("hg19", "hg38")]
    selected = Path(value)
    for name in options:
        bundled = source_root / name
        expected = source_root.resolve() / name
        if (key == "hgnc_table" and name == "data/hgnc/non_alt_loci_set.tsv"
                and bundled.is_symlink() and os.readlink(bundled) == Path(BUNDLED_PATHS[key]).name):
            expected = source_root.resolve() / BUNDLED_PATHS[key]
        legacy_basename = (key in ("clinvar_gene_stat", "interpro_mapping_pickle")
                           and len(selected.parts) == 1 and selected.name == bundled.name
                           and not selected.is_file())
        alias = selected == Path(name) or legacy_basename
        same_source = selected.absolute().resolve() == bundled.resolve()
        if key == "hgnc_table":
            logical = selected if selected.is_absolute() else source_root / selected
            same_source = Path(os.path.abspath(logical)) == Path(os.path.abspath(bundled))
        if ((alias or same_source)
                and bundled.resolve() == expected
                and bundled.is_file() and bundled.stat().st_size > 0):
            return bundled
    return None


def clinvar_default_marker(value: str, source_root: Path | None) -> Path | None:
    """Keep the genuine baseline identity until refresh publishes its candidate.

    Missing/empty defaults are installation inputs to the existing ClinVar
    refresh decision. They are never new output paths inside the checkout.
    """
    if source_root is None:
        return None
    relative = Path("data/ClinVar/clinvar_2star_stats.pkl.gz")
    baseline = source_root / relative
    selected = Path(value)
    legacy = selected == Path(baseline.name) and not selected.is_file()
    if ((selected == relative or selected.absolute() == baseline.absolute() or legacy)
            and baseline.resolve() == source_root.resolve() / relative):
        return baseline
    return None


def spliceai_inspector():
    global _SPLICEAI_INSPECTOR
    if "_SPLICEAI_INSPECTOR" not in globals():
        spec = importlib.util.spec_from_file_location("priva_cache_source_header", Path(__file__).with_name("spliceai_cache_input.py"))
        _SPLICEAI_INSPECTOR = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(_SPLICEAI_INSPECTOR)
    return _SPLICEAI_INSPECTOR


def gnomad_inspector():
    global _GNOMAD_INSPECTOR
    if "_GNOMAD_INSPECTOR" not in globals():
        spec = importlib.util.spec_from_file_location("priva_gnomad_source_family", Path(__file__).with_name("gnomad_cache.py"))
        _GNOMAD_INSPECTOR = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(_GNOMAD_INSPECTOR)
    return _GNOMAD_INSPECTOR


def spliceai_source_format(key: str, value: str) -> str | None:
    if key not in ("spliceai_snv_prescore", "spliceai_indel_prescore") or not Path(value).is_file():
        return None
    # The existing bounded decompressed header inspector detects compression
    # and schema independently of basename/extension. It does not read rows.
    return spliceai_inspector().inspect_header(value)["format"]


def native_spliceai_source(key: str, value: str) -> bool:
    # Its manifest-aware importer owns all 47 shards, including renamed inputs.
    return spliceai_source_format(key, value) == "spliceai2_tsv"


def prepared_spliceai_source(key: str, value: str) -> bool:
    if spliceai_source_format(key, value) != "vcf_spliceai2":
        return False
    if not (Path(value).resolve().parent / "adapter.ready.json").is_file():
        raise ValueError("Prepared SpliceAI2 source requires its completed adjacent adapter.ready.json; ordinary VCF copying is not a substitute")
    return True


def source_is_needed(key: str, config: dict) -> bool:
    if key == "hg38_hg19_chain":
        return config.get("assembly") in ("hg19", "GRCh37")
    if key in ("pext_mean_bw", "pext_tissue_bw_dir"):
        return str(config.get("pext_enabled", "false")).lower() == "true"
    if key == "mavedb_file" and config.get("assembly") in ("hg19", "GRCh37"):
        return False  # The selected GRCh37 profile deliberately omits MaveDB.
    return True


def defer_missing_build_source(key: str, source: Path, source_root: Path | None) -> bool:
    if source_root is None or key not in GOFCARDS_BUILD_ONLY_KEYS:
        return False
    relative = Path(BUNDLED_PATHS[key])
    if source == relative or Path(os.path.abspath(source)) == Path(os.path.abspath(source_root / relative)):
        return False  # Genuine bundled defaults still require existence/nonzero.
    return (key in GOFCARDS_BUILD_ONLY_KEYS and not source.exists()
            and not source.is_symlink()
            and not any(os.environ.get(flag) == "1" for flag in ("PRIVA_FORCE_GOFCARDS_CACHE", "PRIVA_FORCE_ALL_CACHES")))


def automatic_source_paths(config: dict, updates: dict[str, str], source_root: Path) -> dict[str, set[Path]]:
    """Exact public acquisition destinations; arbitrary missing selections fail."""
    if config.get("assembly") not in ("hg19", "GRCh37", "hg38", "GRCh38"):
        return {}
    effective = dict(config, **updates)
    defaults_config = dict(effective)
    for key in resource_groups(config, {key: updates.get(key, config.get(key)) for _, key in FOLDERS.values()}):
        if key not in GENERATED_DIRECTORY_KEYS:
            defaults_config.pop(key, None)
    defaults = _portable_defaults(defaults_config, source_root)
    paths = {key: {Path(defaults[key])} for key in AUTO_SOURCE_KEYS if key in defaults}
    root = checked_storage(config["cache_root"], "cache_root")
    ucsc = "hg19" if config["assembly"] in ("hg19", "GRCh37") else "hg38"
    paths["ref_genome"] = {root / "reference" / f"ucsc.{ucsc}.fasta"}
    # These existing routes preserve the profile's raw assembly spelling.
    # Accept both documented aliases at their standard managed destinations.
    plugin = Path(effective["vep_plugins_cachedir"])
    plugins = Path(effective["vep_plugins_dir"])
    aliases = ("hg19", "GRCh37") if ucsc == "hg19" else ("hg38", "GRCh38")
    for key, name in (("human_ancestor_fasta", "human_ancestor.fa.gz"),
                      ("loftee_conservation_file", "phylocsf_gerp.sql"),
                      ("gerp_bigwig", "gerp_conservation_scores.homo_sapiens.bw")):
        if key != "gerp_bigwig" or ucsc == "hg38":
            paths.setdefault(key, set()).update(plugin / "LoFtee" / alias / name for alias in aliases)
    paths.setdefault("loftee_repo", set()).update(plugins / "LoFtee_repos" / f"loftee-{alias}" for alias in aliases)
    if ucsc == "hg38":
        paths.setdefault("mavedb_file", set()).update(plugin / "MaveDB" / f"MaveDB_variants.{alias}.tsv.gz" for alias in aliases)
    return paths


def resolve_cache_paths(config: dict, source_root: Path | None = None) -> tuple[dict[str, str], dict[str, str]]:
    """Plan destinations from original source paths; never perform adoption here."""
    if not isinstance(config, dict):
        raise ValueError("The configuration must be a YAML mapping")
    validate_repeat_filename(config.get("repeat_region_file_name"))
    if config.get("cache_root") in EMPTY:
        return {}, {}
    root = checked_storage(config["cache_root"], "cache_root")
    folders = config.get("cache_folders")
    if folders is None:
        folders = {}
    if not isinstance(folders, dict) or set(folders) - set(FOLDERS):
        raise ValueError("cache_folders must be a mapping with only the six supported resource names")
    if set(config) & set(FOLDERS):
        raise ValueError("Indent resource folder names beneath cache_folders")
    directories, previous = {}, []
    for resource, (default, key) in FOLDERS.items():
        name = folders.get(resource, default)
        if not isinstance(name, str) or not re.fullmatch(r"[A-Za-z0-9_][A-Za-z0-9_.-]*", name):
            raise ValueError(f"cache_folders.{resource} must be one folder name")
        directories[key] = str(root / name)
        old = config.get(key)
        if isinstance(old, str) and old.startswith("/") and Path(old) != Path("/"):
            previous.append((Path(os.path.normpath(old)), Path(directories[key])))
    previous.sort(key=lambda pair: len(pair[0].parts), reverse=True)
    snv, indel = (config.get(key) for key in ("spliceai_snv_prescore", "spliceai_indel_prescore"))
    if (config.get("assembly") is not None and isinstance(snv, str) and isinstance(indel, str)
            and Path(snv).is_file() and Path(indel).is_file()):
        spliceai_inspector().validate_pair(config["assembly"], snv, indel)
    updates = dict(directories)
    # Mapping tables are selected scientific inputs, not managed cache outputs.
    # Keep explicit cache-side versions and their adjacent receipts in place.
    for key in ("nm_to_enst_tsv", "hgnc_table"):
        mapping = config.get(key)
        if mapping not in EMPTY:
            if not isinstance(mapping, str):
                raise ValueError(f"{key} must be a path string")
            if source_root is not None and not Path(mapping).is_absolute():
                updates[key] = str(Path(os.path.abspath(source_root / mapping)))
    for key, folder in resource_groups(config, directories).items():
        value = config.get(key)
        if value in EMPTY:
            continue
        if not isinstance(value, str):
            raise ValueError(f"{key} must be a path string")
        if key == "hg38_hg19_chain":
            value = str(reverse_chain_source(value, source_root))
        if not source_is_needed(key, config) or defer_missing_build_source(key, Path(value), source_root):
            continue
        if native_spliceai_source(key, value):
            continue
        if prepared_spliceai_source(key, value):
            updates[key] = str(Path(directories["vep_plugins_cachedir"]) / "SpliceAI2/GRCh38" / Path(value).resolve().name)
            continue
        baseline = clinvar_default_marker(value, source_root) if key == "clinvar_gene_stat" else None
        if baseline is not None:
            updates[key] = str(baseline)
            continue
        bundled = bundled_resource(key, value, source_root)
        if bundled is not None:
            updates[key] = str(bundled)
            continue
        path = Path(os.path.abspath(value))
        if key in ("alphamissense_prescore", "gnomad_vcf_chrX") and not path.is_relative_to(root):
            # Keep explicit supplied status through relocation. The resource
            # installer validates these bytes before selecting canonical data;
            # relocation must not turn a bad supplied input into auto-recovery.
            updates[key] = str(folder / "supplied" / key / path.name)
            continue
        # A selected bare output filename stays a generated managed filename.
        for old, new in previous:
            if path.is_relative_to(old):
                updates[key] = str(new / path.relative_to(old))
                break
        else:
            if path.is_relative_to(root):
                continue
            # Names and bytes identify the supplied dataset/release. Keep its
            # basename, with a role-specific directory to avoid collisions.
            updates[key] = str(folder / "supplied" / key / path.name)
    return directories, updates


def adoption_plan(config: dict, updates: dict[str, str], source_root: Path) -> list[dict[str, str]]:
    """Capture explicit source identities before any configuration remapping."""
    for key in ("nm_to_enst_tsv", "hgnc_table"):
        mapping = config.get(key)
        if mapping not in EMPTY:
            if not isinstance(mapping, str):
                raise ValueError(f"{key} must be a path string")
            selected = Path(mapping) if Path(mapping).is_absolute() else source_root / mapping
            if not selected.is_file() or selected.stat().st_size == 0:
                raise ValueError(f"Explicit {key} source is missing/empty; refusing to substitute another dataset/release: {selected}")
    plan = []
    cache_keys = resource_groups(config, {key: updates.get(key, config.get(key)) for _, key in FOLDERS.values()})
    automatic = automatic_source_paths(config, updates, source_root)
    candidates = dict(updates)
    candidates.update({key: updates.get(key, value) for key, value in config.items() if key in cache_keys})
    for key, destination in candidates.items():
        if key not in cache_keys and key not in DIRECTORY_CACHE_KEYS:
            continue  # Source scripts/environment specifications are not caches.
        value = config.get(key)
        if key == "hg38_hg19_chain" and value not in EMPTY:
            value = str(reverse_chain_source(value, source_root))
        if value in EMPTY or not isinstance(value, str) or native_spliceai_source(key, value):
            continue
        if not source_is_needed(key, config):
            continue
        if key == "clinvar_gene_stat" and clinvar_default_marker(value, source_root) is not None:
            continue
        if bundled_resource(key, value, source_root) is not None:
            continue  # No bundle content read/hash/copy, including alias forms.
        source = Path(os.path.abspath(value))
        target = Path(destination)
        if defer_missing_build_source(key, source, source_root):
            continue
        present = source.is_dir() if key in DIRECTORY_CACHE_KEYS else source.is_file()
        if not present:
            if source.exists() or source.is_symlink():
                raise ValueError(f"Explicit {key} source has the wrong file type or is a broken link: {source}")
            if (key in GENERATED_DIRECTORY_KEYS
                    or (key in GENERATED_CACHE_KEYS and key != "clinvar_gene_stat")
                    or source in automatic.get(key, set())):
                continue
            raise ValueError(f"Explicit {key} source is missing; refusing to substitute another dataset/release: {source}")
        if key == "gnomad_vcf_chrX":
            helper = gnomad_inspector()
            canonical = helper.default_path(cache_dir=updates["gnomad_vcf_dir"], assembly=config["assembly"])
            if source != canonical:
                helper.validate_family(chr_x=source, assembly=config["assembly"], reference=config.get("ref_genome"))
            if source != target and source.resolve() != target.resolve():
                original_family, target_family = helper.family_paths(source), helper.family_paths(target)
                for chromosome, member in original_family.items():
                    member_target = target_family[chromosome]
                    plan.append({"key": key, "source": str(member), "destination": str(member_target), "kind": "file"})
                    for suffix in FILE_COMPANIONS:
                        companion = Path(str(member) + suffix)
                        if companion.is_file():
                            plan.append({"key": key, "source": str(companion), "destination": str(member_target) + suffix, "kind": "file"})
            continue
        if source == target or source.resolve() == target.resolve():
            continue
        if prepared_spliceai_source(key, value):
            item = {"key": "spliceai2_prepared_pair", "source": str(source.resolve().parent),
                    "destination": str(target.parent), "kind": "prepared_projection"}
            if item not in plan:
                plan.append(item)
            continue
        if key in DIRECTORY_CACHE_KEYS:
            if not source.is_dir():
                continue  # This is an install destination, not a missing dataset.
            if source.resolve().is_relative_to(source_root.resolve() / "data"):
                continue  # Never recursively hash/copy bundled checkout data.
            plan.append({"key": key, "source": str(source), "destination": str(target), "kind": "directory"})
        elif source.is_file():
            plan.append({"key": key, "source": str(source), "destination": str(target), "kind": "file"})
            for suffix in FILE_COMPANIONS:
                companion = Path(str(source) + suffix)
                if companion.is_file():
                    plan.append({"key": key, "source": str(companion), "destination": str(target) + suffix, "kind": "file"})
        elif key not in GENERATED_CACHE_KEYS or key == "clinvar_gene_stat":
            raise ValueError(f"Explicit {key} source is missing; refusing to substitute another dataset/release: {source}")
    # Parent snapshots already include selected child paths and companions.
    directories = sorted((item for item in plan if item["kind"] == "directory"), key=lambda item: len(Path(item["source"]).parts))
    result = []
    for item in plan:
        if any(item["kind"] != "prepared_projection" and item != parent and Path(item["source"]).is_relative_to(Path(parent["source"]))
               and Path(item["destination"]) == Path(parent["destination"]) / Path(item["source"]).relative_to(Path(parent["source"]))
               for parent in directories):
            continue
        if item not in result:
            result.append(item)
    return sorted(result, key=lambda item: item["kind"] != "prepared_projection")


def cache_identity(path: Path) -> tuple[int, int, int, int]:
    value = path.stat()
    return value.st_ino, value.st_size, value.st_mtime_ns, value.st_ctime_ns


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def adopt_sources(plan: list[dict[str, str]], root: Path) -> None:
    """Atomic byte-preserving copies; never mutate a source or use hardlinks."""
    if not plan:
        return
    root.mkdir(parents=True, exist_ok=True)
    records = []
    prepared_sources = {Path(item["source"]).resolve() for item in plan if item["kind"] == "prepared_projection"}
    with (root / ".installer-adoption.lock").open("a") as lock:
        fcntl.flock(lock, fcntl.LOCK_EX)

        def require_managed_parent(destination):
            ancestor = destination.parent
            while not ancestor.exists() and not ancestor.is_symlink():
                ancestor = ancestor.parent
            if not ancestor.resolve().is_relative_to(root.resolve()):
                raise ValueError(f"Refusing to copy through a managed-directory link into an external source: {destination}")

        def copy_file(source, destination):
            source, destination = Path(source), Path(destination)
            before = cache_identity(source)
            require_managed_parent(destination)
            destination.parent.mkdir(parents=True, exist_ok=True)
            descriptor, name = tempfile.mkstemp(prefix=".adopting-", dir=destination.parent)
            temporary = Path(name)
            os.close(descriptor)
            try:
                shutil.copy2(source, temporary)
                source_sha = sha256_file(source)
                if cache_identity(source) != before or sha256_file(temporary) != source_sha:
                    raise ValueError(f"Source changed or copy differs during cache adoption: {source}")
                if destination.exists() and not destination.is_symlink() and destination.is_file() and sha256_file(destination) == source_sha:
                    temporary.unlink()
                else:
                    if destination.exists() or destination.is_symlink():
                        stamp = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%S.%fZ")
                        retained = root / ".installer-retained" / stamp / destination.relative_to(root)
                        retained.parent.mkdir(parents=True, exist_ok=True)
                        destination.rename(retained)
                    temporary.replace(destination)
                records.append({"source": str(source), "destination": str(destination), "bytes": before[1],
                                "sha256": source_sha, "source_identity": list(before)})
                return str(destination)
            finally:
                temporary.unlink(missing_ok=True)

        for item in plan:
            source, destination = Path(item["source"]), Path(item["destination"])
            if not destination.is_relative_to(root):
                raise ValueError(f"Adoption destination escapes cache_root: {destination}")
            if item["kind"] == "prepared_projection":
                # The existing SpliceAI2 receipt validator/relocator owns this
                # physical copy. Generic copying cannot preserve absolute-path
                # receipt contracts or prove projection reuse.
                require_managed_parent(destination)
                spliceai_inspector().adopt_completed_projection(
                    source_dir=source, cache_dir=destination, cache_root=root)
                records.append({"source": str(source), "destination": str(destination),
                                "kind": "prepared_projection", "source_identity": list(cache_identity(source))})
            elif item["kind"] == "directory":
                if destination.resolve().is_relative_to(source.resolve()):
                    raise ValueError(f"Cannot adopt a cache directory into itself: {source}")
                require_managed_parent(destination)
                if destination.is_symlink():
                    raise ValueError(f"Refusing mutable directory adoption through an existing link: {destination}")
                if destination.is_dir():
                    for parent, directories, _ in os.walk(destination, followlinks=False):
                        for name in directories:
                            child = Path(parent) / name
                            if child.is_symlink():
                                raise ValueError(f"Refusing mutable cache adoption through an existing directory link: {child}")
                # Snapshot copies preserve originals and all directory-local
                # indexes/receipts. Validate every copied source's stat identity
                # again before publishing the configuration or adoption record.
                def preserve_completed_projection(directory, names):
                    return [name for name in names if (Path(directory) / name).resolve() in prepared_sources]
                shutil.copytree(source, destination, symlinks=False, dirs_exist_ok=True,
                                copy_function=copy_file, ignore=preserve_completed_projection)
            else:
                copy_file(source, destination)
        if any(list(cache_identity(Path(item["source"]))) != item["source_identity"] for item in records):
            raise ValueError("An explicit cache source changed during adoption; original config is retained")
        receipt = root / "installer-source-adoptions.json"
        previous = json.loads(receipt.read_text()) if receipt.is_file() else {
            "schema_version": 1, "copies": [], "scientific_acceptance": False,
            "scope": "Byte-preserving source adoption; resource-specific scientific reuse checks remain required."}
        previous["copies"].extend(record for record in records if record not in previous["copies"])
        descriptor, temporary = tempfile.mkstemp(prefix=".adoptions-", dir=root)
        with os.fdopen(descriptor, "w") as handle:
            json.dump(previous, handle, indent=2)
            handle.write("\n")
        os.replace(temporary, receipt)


class UniqueLoader(yaml.SafeLoader):
    """Do not let duplicate YAML keys silently change a selected input."""


def unique_mapping(loader: UniqueLoader, node: yaml.MappingNode, deep: bool = False) -> dict:
    result = {}
    for key_node, value_node in node.value:
        key = loader.construct_object(key_node, deep=deep)
        if key in result:
            raise ValueError(f"Duplicate configuration key: {key}")
        result[key] = loader.construct_object(value_node, deep=deep)
    return result


UniqueLoader.add_constructor(yaml.resolver.BaseResolver.DEFAULT_MAPPING_TAG, unique_mapping)


def checked_storage(value: object, key: str) -> Path:
    if not isinstance(value, str) or not re.fullmatch(r"/[A-Za-z0-9_./-]+", value):
        raise ValueError(f"Set {key} to an absolute storage path without spaces or shell characters")
    path = Path(os.path.normpath(value))
    if path == Path("/"):
        raise ValueError(f"{key} must name a dedicated storage directory")
    return path


def bundled_cache_default(key: str, fallback: Path, source_root: Path, *, relative: str | None = None) -> Path:
    """Known bundled defaults are read-only; force/new outputs use managed storage."""
    relative = relative or BUNDLED_CACHES[key]
    bundled = source_root / relative
    forced = any(os.environ.get(flag) == "1" for flag in ("PRIVA_FORCE_ALL_CACHES", *FORCE_KEYS.get(key, ())))
    if (not forced and bundled.resolve() == source_root.resolve() / relative
            and bundled.is_file() and bundled.stat().st_size):
        return bundled
    return fallback


def _portable_defaults(config: dict, source_root: Path) -> dict[str, str]:
    """Fill unset fields after selected paths have a managed destination."""
    if not isinstance(config, dict):
        raise ValueError("The configuration must be a YAML mapping")
    root = checked_storage(config.get("cache_root"), "cache_root")
    checked_storage(config.get("tmp_dir"), "tmp_dir")
    assembly = {"hg19": "hg19", "GRCh37": "hg19", "hg38": "hg38", "GRCh38": "hg38"}.get(config.get("assembly"))
    if assembly is None:
        raise ValueError("Set assembly to hg19/GRCh37 or hg38/GRCh38")
    if not re.fullmatch(r"[1-9][0-9]*", str(config.get("vep_installed_version", ""))):
        raise ValueError("Set vep_installed_version to the requested numeric major release")
    folders = config.get("cache_folders") or {}
    if not isinstance(folders, dict) or set(folders) - set(FOLDERS):
        raise ValueError("cache_folders must contain only the six supported resource names")
    updates = {}
    effective = dict(config)

    def fill(key: str, value: Path | str) -> None:
        if effective.get(key) in EMPTY:
            updates[key] = str(value)
            effective[key] = str(value)

    for resource, (name, key) in FOLDERS.items():
        selected_name = folders.get(resource, name)
        if not isinstance(selected_name, str) or not re.fullmatch(r"[A-Za-z0-9_][A-Za-z0-9_.-]*", selected_name):
            raise ValueError(f"cache_folders.{resource} must be one folder name")
        fill(key, root / selected_name)
    for key, relative in BUNDLED_PATHS.items():
        selected = effective.get(key)
        if selected in EMPTY:
            fill(key, source_root / relative)
        elif not isinstance(selected, str):
            raise ValueError(f"{key} must be a path string")
        elif not Path(selected).is_absolute():
            updates[key] = str(source_root / selected)
            effective[key] = updates[key]

    plugin = Path(effective["vep_plugins_cachedir"])
    clinvar = Path(effective["clinvar_vcf_dir"])
    genome = "GRCh37" if assembly == "hg19" else "GRCh38"
    # LoFtee's existing installer uses the configured spelling for its directories.
    loftee_assembly = config["assembly"]
    alpha = plugin / "AlphaMissense"
    vcf = clinvar / f"clinvar.{assembly}.vep.vcf.gz"
    fill("gnomad_vcf_chrX", Path(effective["gnomad_vcf_dir"]) / (
        "gnomad.joint.v4.1.sites.hg19.chrX.vcf.gz" if assembly == "hg19"
        else "gnomad.joint.v4.1.sites.chrX.vcf.bgz"))
    fill("clinvar_vcf", vcf)
    for key, suffix in (("clinvar_aa_stat", "aa_change.stats.pkl"),
                        ("clinvar_splice_stat", "splice_change.stats.pkl"),
                        ("clinvar_patho_af_stat", "patho_af_stat.pkl"),
                        ("clinvar_patho_exon_af_stat", "patho_exon_af_stat.pkl")):
        fill(key, str(effective["clinvar_vcf"]).removesuffix(".vcf.gz") + "." + suffix)
    # Preserve this baseline marker even when absent/empty. The existing
    # ClinVar refresh procedure recognizes it and selects a managed candidate.
    fill("clinvar_gene_stat", source_root / "data/ClinVar/clinvar_2star_stats.pkl.gz")
    fill("utr_annotator_file", plugin / f"UTRannotator/uORF_5UTR_{genome}_PUBLIC.txt")
    fill("loeuf_prescore", source_root / "data/loeuf" / (
        "loeuf_dataset.tsv.gz" if assembly == "hg19" else "loeuf_dataset_hg38.tsv.gz"))
    fill("loftee_parent_dir", Path(effective["vep_plugins_dir"]) / "LoFtee_repos")
    fill("loftee_repo", Path(effective["loftee_parent_dir"]) / f"loftee-{loftee_assembly}")
    loftee = plugin / "LoFtee" / loftee_assembly
    fill("human_ancestor_fasta", loftee / "human_ancestor.fa.gz")
    fill("loftee_conservation_file", loftee / "phylocsf_gerp.sql")
    if assembly == "hg38":
        fill("gerp_bigwig", loftee / "gerp_conservation_scores.homo_sapiens.bw")
    for key, name in {
        "alphamissense_prescore": f"AlphaMissense_{assembly}.tsv.gz",
        "alphamissense_vcf": f"AlphaMissense_{assembly}.vcf.gz",
        "alphamissense_vep_vcf": f"AlphaMissense_{assembly}.vep.vcf.gz",
        "alphamissense_vep_base_vcf": f"AlphaMissense_{assembly}.vep.base.vcf.gz",
        "alphamissense_pd_stat": f"AlphaMissense_{assembly}.vep.prot.domain.stats.pkl",
        "alphamissense_tranx_domain_map": f"transcript_pp_domain_mapping.{assembly}.pkl",
        "alphamissense_intolerant_domains": f"domain_tolerance_analysis.{assembly}.tsv",
        "alphamissense_intolerant_motifs": f"AlphaMissense_{assembly}.kde.pkl",
        "all_intolerant_domains": f"all_intolerant_domains.{assembly}.pkl",
        "pm1_regions_pkl": f"pm1_regions.{assembly}.pkl",
    }.items():
        fill(key, alpha / name)
    fill("conservation_file", plugin / f"Conservation/gerp_conservation_scores.homo_sapiens.{genome}.bw")
    fill("splicevault_prescore", plugin / f"SpliceVault/SpliceVault_data_{genome}.tsv.gz")
    if assembly == "hg38":
        fill("mavedb_file", plugin / f"MaveDB/MaveDB_variants.{config['assembly']}.tsv.gz")
    if str(config.get("pext_enabled", "false")).lower() == "true":
        fill("pext_mean_bw", plugin / "pext" / assembly / ("mean_proportion.bw" if assembly == "hg19" else "exp_prop_mean.bw"))
        fill("pext_tissue_bw_dir", plugin / "pext" / assembly)
    fill("cadd_base_dir", Path(effective["cadd_parent_dir"]) / "CADD-scripts-1.7.2")
    fill("hg38_hg19_chain", source_root / "data/liftover/hg38ToHg19.over.chain.gz")
    fill("repeat_region_file_name", f"ucsc.{assembly}.fasta.repeat.bed.gz")
    bundled_interpro = source_root / "data/InterPro/Interpro_entry_mapping.pkl.gz"
    # This known bundled table is already supplied by the checkout. Do not
    # manufacture a replacement download or inspect external symlink contents.
    if (bundled_interpro.resolve().is_relative_to(source_root.resolve() / "data/InterPro")
            and bundled_interpro.is_file() and bundled_interpro.stat().st_size):
        fill("interpro_mapping_pickle", bundled_interpro)
    else:
        fill("interpro_mapping_pickle", root / "InterPro/Interpro_entry_mapping.pkl.gz")
    fill("clingen_map", bundled_cache_default("clingen_map", root / f"ClinGen/clingen_map.{assembly}.pkl.gz", source_root,
                                               relative=f"data/clingen/clingen_map.{assembly}.pkl.gz"))
    fill("gene_dosage_sensitivity", bundled_cache_default("gene_dosage_sensitivity", root / f"ClinGen/gene_dosage_sensitivity.{assembly}.tsv", source_root,
                                                           relative=f"data/clingen/gene_dosage_sensitivity.{assembly}.tsv"))
    fill("cds_fasta_file", plugin / f"Ensembl_CDS/release-{config['vep_installed_version']}/Homo_sapiens.{genome}.cds.all.fa.gz")

    fill("gene_mechanism_cache_dir", root / "patho_mechanism")
    mechanism = Path(effective["gene_mechanism_cache_dir"])
    fill("gene_mechanism_raw_dir", mechanism / "raw")
    fill("clinvar_vcv_match_audit", mechanism / "clinvar_vcv_gofcards_matches.tsv")
    fill("gene_nonlof_mechanism_json", bundled_cache_default("gene_nonlof_mechanism_json", mechanism / "gene_nonlof_mechanism_curated_assertions.json.gz", source_root))
    fill("hpo_condition_mechanism_json", bundled_cache_default("hpo_condition_mechanism_json", mechanism / "hpo_condition_mechanism_cache.json.gz", source_root))
    fill("ddg2p_mechanism_evidence", bundled_cache_default("ddg2p_mechanism_evidence", mechanism / "gene_pathogenic_mechanism_evidence.tsv", source_root))
    fill("gofcards_exact_gof_cache", bundled_cache_default("gofcards_exact_gof_cache", root / "gofcards/gofcards_exact_gof.json.gz", source_root))
    fill("gofcards_workdir", root / "gofcards/build_work")
    fill("clinvar_vcv_xml", Path(effective["gene_mechanism_raw_dir"]) / "clinvar_vcv/ClinVarVCVRelease_00-latest_weekly.xml.gz")
    hpo_release = config.get("hpo_release") or "v2026-06-23"
    mondo_release = config.get("mondo_release") or "v2026-07-06"
    fill("hpo_raw_dir", root / "hpo/raw" / hpo_release)
    fill("hpo_assertions", root / "hpo/genes_to_phenotype.assertions.tsv.gz")
    fill("mondo_obo", root / f"mondo/mondo-simple.{mondo_release}.obo")
    fill("mondo_disease_scope_registry", root / "mondo/disease_scope.tsv.gz")
    # Reference/chain and NM mapping planners retain their own source contracts.
    # Licensed inputs, MaveDB/pext and workflow-only paths remain unset here.
    return updates


def portable_updates(config: dict, source_root: Path) -> dict[str, str]:
    """Plan from the original selections; acquisition/adoption is separate."""
    _, updates = resolve_cache_paths(config, source_root)
    effective = dict(config, **updates)
    updates.update(_portable_defaults(effective, source_root))
    # Fill-generated directory defaults must also be managed when legacy config
    # explicitly selected an old output/raw directory outside the cache root.
    return updates


def render_config(text: str, updates: dict[str, str]) -> str:
    """Keep comments and all unrelated scalars intact."""
    if text.lstrip().startswith("{"):
        # JSON is also a YAML mapping; keep it JSON for direct helper callers.
        config = json.loads(text)
        config.update(updates)
        return json.dumps(config, indent=2, ensure_ascii=False) + "\n"
    for key, value in sorted(updates.items()):
        pattern = re.compile(r"^" + re.escape(key) + r"[ \t]*:[^\n]*$", re.MULTILINE)
        matches = list(pattern.finditer(text))
        if len(matches) > 1:
            raise ValueError(f"Duplicate configuration key: {key}")
        comment = ""
        if matches:
            old = matches[0].group().split(":", 1)[1].strip()
            scalar = re.fullmatch(r'''(?:"(?:\\.|[^"\\])*"|'(?:''|[^'])*'|[^#]*?)(\s+#.*)?''', old)
            if scalar and scalar.group(1):
                comment = scalar.group(1)
        replacement = key + ": " + json.dumps(value) + comment
        text = pattern.sub(lambda _: replacement, text) if matches else text.rstrip("\n") + "\n" + replacement + "\n"
    return text


def resolve_config(args: argparse.Namespace) -> None:
    original = args.config.read_text()
    config = yaml.load(original, Loader=UniqueLoader)
    managed = json.loads(args.managed_updates_json) if args.managed_updates_json else {}
    if not isinstance(managed, dict) or not all(isinstance(key, str) and isinstance(value, str)
                                               for key, value in managed.items()):
        raise ValueError("Managed path updates must be a JSON mapping of path strings")
    source_root = args.source_root.absolute()
    updates = portable_updates(config, source_root)
    if any(key in updates and updates[key] != value for key, value in managed.items()):
        raise ValueError("Managed updates differ from the original configuration/cache_root plan")
    plan = adoption_plan(config, updates, source_root)
    if args.check:
        print(json.dumps(updates, indent=2, sort_keys=True))
        return
    if args.output is not None:
        # Prerequisites must see the selected source identity, including a
        # genuinely missing custom source. No bytes are adopted during planning.
        for key in resource_groups(config, {key: updates[key] for _, key in FOLDERS.values()}):
            value = config.get(key)
            if key == "hg38_hg19_chain" and value not in EMPTY:
                value = str(reverse_chain_source(value, source_root))
            if isinstance(value, str) and value and key not in {key for _, key in FOLDERS.values()}:
                source = Path(os.path.abspath(value))
                if source.is_file() or source.is_dir() or (key not in GENERATED_CACHE_KEYS and bundled_resource(key, value, source_root) is None):
                    updates[key] = value
    else:
        adopt_sources(plan, checked_storage(config["cache_root"], "cache_root"))
    updated = render_config(original, updates)
    destination = args.output or args.config
    if destination == args.config and updated == original:
        return
    # Source cache bytes/companions have been adopted before this atomic config
    # selection. Scientific reuse/rebuild checks still belong to each installer.
    descriptor, temporary = tempfile.mkstemp(prefix="." + destination.name + ".", dir=destination.parent)
    try:
        with os.fdopen(descriptor, "w") as handle:
            handle.write(updated)
        os.chmod(temporary, args.config.stat().st_mode & 0o777)
        if args.config.read_text() != original:
            raise ValueError("Configuration changed during path resolution; retry after its writer closes")
        os.replace(temporary, destination)
    finally:
        if os.path.exists(temporary):
            os.unlink(temporary)


def configure_cache_paths(config_file, config, mode, source_root=None):
    """Compatibility procedure for the shell's managed-root phase."""
    if mode not in ("write", "--check"):
        raise ValueError("Use priva_configure_cache_paths CONFIG [--check]")
    source_root = Path(source_root) if source_root is not None else Path(__file__).resolve().parents[1]
    directories, updates = resolve_cache_paths(config, source_root)
    if mode == "--check":
        print(json.dumps(updates, indent=2, sort_keys=True))
        return
    if not updates:
        return
    path = Path(config_file)
    original = path.read_text()
    updated = render_config(original, updates)
    adopt_sources(adoption_plan(config, updates, source_root), checked_storage(config["cache_root"], "cache_root"))
    for directory in directories.values():
        Path(directory).mkdir(parents=True, exist_ok=True)
    if updated != original:
        descriptor, temporary = tempfile.mkstemp(prefix="." + path.name + ".", dir=path.parent)
        try:
            with os.fdopen(descriptor, "w") as stream:
                stream.write(updated)
            os.chmod(temporary, path.stat().st_mode & 0o777)
            if path.read_text() != original:
                raise ValueError("Configuration changed during source adoption; original selections retained")
            os.replace(temporary, path)
        finally:
            if os.path.exists(temporary):
                os.unlink(temporary)
    for key, directory in directories.items():
        print(f"{key}: {directory}")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    commands = parser.add_subparsers(dest="command", required=True)
    resolve = commands.add_parser("resolve", help="Fill unset managed paths and resolve checkout-relative bundled paths")
    resolve.add_argument("--config", type=Path, required=True)
    resolve.add_argument("--source-root", type=Path, required=True)
    resolve.add_argument("--check", action="store_true", help="Print planned updates without changing files or directories")
    resolve.add_argument("--managed-updates-json", help="Existing shell managed-root resolver's JSON updates")
    resolve.add_argument("--output", type=Path, help="Write an effective preflight config here, retaining the selected config")
    managed = commands.add_parser("managed-roots", help="Adopt selected cache sources and configure managed resource directories")
    managed.add_argument("--config", type=Path, required=True)
    managed.add_argument("--source-root", type=Path, required=True)
    managed.add_argument("--check", action="store_true", help="Print destination updates without copying or writing")
    return parser.parse_args()


def main() -> int:
    try:
        args = parse_args()
        if args.command == "managed-roots":
            config = yaml.load(args.config.read_text(), Loader=UniqueLoader)
            configure_cache_paths(args.config, config, "--check" if args.check else "write", args.source_root)
        else:
            resolve_config(args)
    except (OSError, ValueError, TypeError, yaml.YAMLError) as error:
        print(f"ERROR: {error}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
