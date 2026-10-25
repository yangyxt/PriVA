#!/usr/bin/env python3
"""Installer input checks; manual-input prompts can atomically fill unset paths."""
from __future__ import annotations

import argparse
import csv
import gzip
import importlib.util
import os
from pathlib import Path
import pickle
import re
import sys
import tempfile

import yaml

_MAPPING_SPEC = importlib.util.spec_from_file_location(
    "priva_biomart_mapping", Path(__file__).with_name("biomart_nm_to_enst.py"))
BIOMART = importlib.util.module_from_spec(_MAPPING_SPEC)
_MAPPING_SPEC.loader.exec_module(BIOMART)
_REFERENCE_SPEC = importlib.util.spec_from_file_location(
    "priva_reference_resources", Path(__file__).with_name("reference_resources.py"))
REFERENCES = importlib.util.module_from_spec(_REFERENCE_SPEC)
_REFERENCE_SPEC.loader.exec_module(REFERENCES)
_CONFIG_SPEC = importlib.util.spec_from_file_location(
    "priva_installer_config_paths", Path(__file__).with_name("installer_config_paths.py"))
CONFIG_PATHS = importlib.util.module_from_spec(_CONFIG_SPEC)
_CONFIG_SPEC.loader.exec_module(CONFIG_PATHS)


WORKFLOW_ONLY_KEYS = {
    "input_vcf", "ped_file", "output_dir", "hub_file_dir", "alt_disease_vcf",
    "control_vcf", "pp1_vcf", "pp1_ped", "conda_env_yaml",
}
CHR1_LENGTH = {"hg19": 249250621, "hg38": 248956422}
MANUAL_INPUT_KEYS = ("primateai_prescore", "spliceai_snv_prescore", "spliceai_indel_prescore")


def checked_path(value: str, key: str, _depth: int = 0) -> Path:
    if not isinstance(value, str) or not value.strip() or value == "null":
        raise ValueError(f"Set {key} to the original, explicitly supplied input path")
    if value.startswith("/paedyl01/"):
        raise ValueError(f"Unresolved old-host path {key}: {value}; supply its actual local path")
    path = Path(os.path.abspath(value))
    # Optional deployment policy: reject paths before inspecting them. The
    # installer itself remains usable on hosts with other filesystem layouts.
    roots = [Path(os.path.abspath(root)) for root in
             os.environ.get("PRIVA_INSTALL_INSPECT_ROOTS", "").split(os.pathsep) if root]
    if roots:
        if _depth > 40:
            raise ValueError(f"Too many symlinks for {key}: {value}")
        matches = [root for root in roots if path == root or root in path.parents]
        if not matches:
            raise ValueError(f"Path outside permitted inspection roots for {key}: {value}")
        root = max(matches, key=lambda item: len(item.parts))
        cursor = root
        for part in path.relative_to(root).parts:
            cursor /= part
            if cursor.is_symlink():
                # Do not follow a link into data outside the selected policy.
                target = Path(os.path.abspath(cursor.parent / os.readlink(cursor)))
                if not any(target == allowed or allowed in target.parents for allowed in roots):
                    raise ValueError(f"Symlink leaves permitted inspection roots for {key}: {value}")
                cursor = checked_path(str(target), key, _depth + 1)
        path = cursor
    return path


def readable_file(value: str, key: str) -> Path:
    path = checked_path(value, key)
    if not path.is_file() or not os.access(path, os.R_OK) or path.stat().st_size == 0:
        raise ValueError(f"Missing, empty, or unreadable {key}: {value}")
    return path


def _unset_manual_input(value: object) -> bool:
    return value is None or isinstance(value, str) and value.strip() in ("", "null")


def validate_manual_file(value: str, key: str) -> None:
    """Check readiness without reading scientific payloads or creating indexes."""
    path = readable_file(value, key)

    def unfinished(part: str) -> bool:
        part = part.lower()
        return (".partial" in part or part == ".building" or part.startswith(".building-")
                or part.startswith(".download.")
                or part.endswith((".part", ".crdownload", ".download", ".copying", ".tmp")))

    # A final-looking symlink must not conceal an explicitly unfinished source.
    if any(unfinished(part) for selected in (Path(os.path.abspath(value)), path.resolve())
           for part in selected.parts):
        raise ValueError(f"Incomplete/staging manual input {key}: {value}; select the completed original file")


def manual_inputs(config_file: Path, *, prompt: bool = False) -> None:
    """Validate all manual selections before acquisition; optionally fill only unset keys.

    This metadata-only gate does not establish assembly, headers, index validity,
    or scientific suitability. Those checks stay with each resource installer.
    Explicit invalid selections are never replaced, including in prompt mode.
    """
    destination = config_file.resolve(strict=True)
    original = destination.read_text(encoding="utf-8")
    config = yaml.load(original, Loader=CONFIG_PATHS.UniqueLoader)
    if not isinstance(config, dict):
        raise ValueError("The configuration must be a YAML mapping")
    missing = []
    for key in MANUAL_INPUT_KEYS:
        if _unset_manual_input(config.get(key)):
            missing.append(key)
        else:
            validate_manual_file(config[key], key)
    if not missing:
        return  # Prepared inputs never cause a config write, even with --prompt.
    if not prompt or not sys.stdin.isatty():
        raise ValueError("Manual publisher inputs must be prepared before resource installation: "
                         + ", ".join(missing)
                         + "; set these paths in config.yaml, or use --stage manual-inputs --prompt "
                           "from an interactive terminal. No files were acquired or configuration changed")
    updates = {}
    print("Select completed manual publisher files; all paths are checked before the configuration is saved.",
          file=sys.stderr, flush=True)
    for key in missing:
        try:
            updates[key] = input(f"{key}: ").strip()
        except EOFError as error:
            raise ValueError("Manual input collection ended; configuration is unchanged") from error
    # Collect all missing paths, then validate them as one transaction. A later
    # invalid answer must not leave earlier answers saved in the configuration.
    selected = dict(config, **updates)
    for key in MANUAL_INPUT_KEYS:
        validate_manual_file(selected[key], key)
    updated = CONFIG_PATHS.render_config(original, updates)
    # The shared renderer handles the distributed template's scalar keys. Do
    # not publish an ambiguous edit if a hand-edited YAML construct cannot be
    # preserved by it (for example an anchor or a quoted duplicate key).
    if yaml.load(updated, Loader=CONFIG_PATHS.UniqueLoader) != selected:
        raise ValueError("Manual input paths could not be updated without changing other configuration values")
    descriptor, temporary = tempfile.mkstemp(prefix="." + destination.name + ".", dir=destination.parent)
    try:
        with os.fdopen(descriptor, "w", encoding="utf-8") as handle:
            handle.write(updated)
            handle.flush()
            os.fsync(handle.fileno())
        os.chmod(temporary, destination.stat().st_mode & 0o777)
        if destination.read_text(encoding="utf-8") != original:
            raise ValueError("Configuration changed during manual input collection; retry after its writer closes")
        os.replace(temporary, destination)
    finally:
        if os.path.exists(temporary):
            os.unlink(temporary)


def validate_nm_to_enst(value: str, assembly: str | None = None, source_root: Path | None = None) -> None:
    readable_file(value, "nm_to_enst_tsv")  # Check policy before opening chosen path or its sidecar.
    path = Path(os.path.abspath(value))
    named = BIOMART.source_for_basename(path.name) if path.name in BIOMART.PINNED_BASENAMES else None
    if named is not None and assembly is not None and assembly not in named["aliases"]:
        raise ValueError(f"Pinned BioMart mapping {path.name} is not compatible with {assembly}")
    if BIOMART.is_bundled_mapping(path, source_root):
        # Existing genuine bundled data follow the user's existence/nonzero
        # policy. External symlink targets and new FTP resources remain checked.
        return
    metadata = checked_path(str(BIOMART.provenance_path(path)), "nm_to_enst_tsv provenance")
    if metadata.exists() or metadata.is_symlink() or path.name in BIOMART.PINNED_BASENAMES:
        BIOMART.verify_pinned(path, assembly)
        return
    # Both unchanged scientific readers use columns 2 and 5. Validate those
    # original columns; never manufacture a replacement mapping or rewrite it.
    with path.open(newline="", encoding="utf-8-sig") as handle:
        rows = csv.reader(handle, delimiter="\t")
        header = next(rows, [])
        if len(header) < 5 or header[1].strip() != "Transcript stable ID" or header[4].strip() != "RefSeq mRNA ID":
            raise ValueError(f"nm_to_enst_tsv must retain the original five-column Ensembl/RefSeq layout: {value}")
        if not any(len(row) >= 5 and row[1].startswith("ENST") and row[4].startswith("NM_") for row in rows):
            raise ValueError(f"nm_to_enst_tsv contains no usable original NM-to-ENST rows: {value}")


def _mapping_default_source(config: dict) -> dict:
    release = config.get("vep_installed_version")
    if release in (None, "", "null"):
        release = None
    return BIOMART.source_for_assembly(config.get("assembly"), release)


def _mapping_default_path(config: dict, source_root: Path, selected: dict) -> Path:
    bundled = Path(os.path.abspath(source_root / selected["relative_path"]))
    if os.path.lexists(bundled):
        return bundled  # Reuse an existing selected source; never overwrite it.
    root = config.get("cache_root")
    if root in (None, "", "null"):
        return Path(os.path.abspath(source_root / selected["relative_path"]))
    if not isinstance(root, str):
        raise ValueError("cache_root must be a path string")
    root = Path(root)
    if not root.is_absolute():
        root = source_root / root
    return Path(os.path.abspath(root / selected["relative_path"]))


def nm_to_enst_path(config: dict, source_root: Path) -> Path:
    """Preserve explicit inputs; plan new defaults from assembly and selected VEP."""
    value = config.get("nm_to_enst_tsv")
    if value in (None, "", "null"):
        selected = _mapping_default_source(config)
        value = str(_mapping_default_path(config, source_root, selected))
    elif not isinstance(value, str):
        raise ValueError("nm_to_enst_tsv must be a path string")
    elif not Path(value).is_absolute():
        value = str(source_root / value)
    # Preserve the logical filename/adjacent provenance when DataLad resolves
    # the payload into .git/annex. checked_path still enforces every link guard.
    checked_path(value, "nm_to_enst_tsv")
    path = Path(os.path.abspath(value))
    named_source = BIOMART.source_for_basename(path.name) if path.name in BIOMART.PINNED_BASENAMES else None
    if named_source is not None and config.get("assembly") not in named_source["aliases"]:
        raise ValueError(
            f"Pinned BioMart mapping {path.name} is not compatible with {config.get('assembly')}; "
            "supply an explicitly assembly-compatible nm_to_enst_tsv"
        )
    return path


def mapping_fetch_source(config: dict, source_root: Path) -> str:
    """Authorize acquisition only for the missing automatically planned destination."""
    target = nm_to_enst_path(config, source_root)
    if target.is_symlink() and not target.exists():
        raise ValueError(f"Explicit nm_to_enst_tsv is missing: {target}; restore that selected linked input")
    if config.get("nm_to_enst_tsv") not in (None, "", "null"):
        raise ValueError(f"Explicit nm_to_enst_tsv is missing: {target}; restore that selected input")
    selected = _mapping_default_source(config)
    bundled = _mapping_default_path(config, source_root, selected)
    if target != bundled:
        raise ValueError(f"Explicit nm_to_enst_tsv is missing: {target}; restore that selected input")
    if os.path.lexists(target):
        raise ValueError(f"Bundled nm_to_enst_tsv already exists; validate or recover it explicitly: {target}")
    metadata = checked_path(str(BIOMART.provenance_path(target)), "nm_to_enst_tsv provenance")
    if os.path.lexists(metadata):
        raise ValueError(f"Missing bundled nm_to_enst_tsv has existing provenance: {metadata}; retain it and recover explicitly")
    return selected["source_id"]


def validate_mapping(config: dict, source_root: Path) -> Path:
    path = nm_to_enst_path(config, source_root)
    validate_nm_to_enst(str(path), config.get("assembly"), source_root)
    return path


class _StatsUnpickler(pickle.Unpickler):
    def find_class(self, module, name):
        if module == "builtins" and name in {"set", "frozenset"}:
            return getattr(__import__(module), name)
        raise ValueError(f"Unexpected object in bundled ClinVar statistics: {module}.{name}")


def validate_clinvar_stats(path: Path) -> None:
    fields = {"ensembl_id", "clinvar_pathogenic", "large_aachange", "small_aachange", "putative_nmd_variants"}
    try:
        readable_file(str(path), "bundled clinvar_gene_stat")
        with gzip.open(path, "rb") as handle:
            stats = _StatsUnpickler(handle).load()
        if not isinstance(stats, dict) or not stats or not all(
                isinstance(value, dict) and fields <= value.keys() for value in stats.values()):
            raise ValueError("missing gene statistics fields")
    except (OSError, ValueError, EOFError, pickle.UnpicklingError) as error:
        raise ValueError(
            f"Cannot reuse bundled ClinVar gene statistics {path}: {error}. "
            "Restore a validated original artifact for baseline reuse. Separate candidate regeneration "
            "requires explicit ClinVar variant_summary and MANE inputs: run "
            "ClinVar_Gene_stat_regenerate CONFIG SUMMARY MANE OUTPUT with a distinct output; "
            "this does not replace the baseline."
        ) from error


def clinvar_stats_path(config: dict, source_root: Path) -> Path:
    selected = str(config.get("clinvar_gene_stat") or "clinvar_2star_stats.pkl.gz")
    path = Path(selected)
    # A bare legacy name is bundled-relative. Retain explicit absolute paths
    # and relative paths with directories, as the runtime consumer does.
    return path if path.is_absolute() or "/" in selected else source_root / "data/ClinVar" / path



def validate_selected_clinvar_stats(config: dict, source_root: Path) -> None:
    """Bundled original data need only exist; supplied/external bytes keep their checks."""
    path = clinvar_stats_path(config, source_root)
    bundled = source_root / "data/ClinVar/clinvar_2star_stats.pkl.gz"
    source_data = source_root.resolve() / "data/ClinVar"
    if path == bundled and path.resolve().is_relative_to(source_data):
        try:
            readable_file(str(path), "bundled clinvar_gene_stat")
        except ValueError as error:
            raise ValueError(
                f"{error}. Restore the nonempty original bundled artifact for baseline reuse. "
                "Separate candidate regeneration requires explicit ClinVar variant_summary and MANE inputs: "
                "run ClinVar_Gene_stat_regenerate CONFIG SUMMARY MANE OUTPUT with a distinct output; "
                "this does not replace the baseline."
            ) from error
    else:
        validate_clinvar_stats(path)


def validate_main(config: dict, source_root: Path, *, allow_missing_bundled_mapping: bool = False,
                  allow_missing_bundled_clinvar: bool = False) -> None:
    problems = [f"Unresolved old-host path {key}: {value}; supply its actual local path"
                for key, value in config.items()
                if key not in WORKFLOW_ONLY_KEYS and isinstance(value, str) and value.startswith("/paedyl01/")]
    # Reject old-host values before opening any of their paths.
    if problems:
        raise ValueError("\n".join(problems))
    mapping = nm_to_enst_path(config, source_root)
    if allow_missing_bundled_mapping and config.get("nm_to_enst_tsv") in (None, "", "null") and not mapping.exists():
        mapping_fetch_source(config, source_root)  # Reject orphan provenance before deferring acquisition.
    else:
        validate_mapping(config, source_root)
    expected = source_root / "data/hpo/genes_to_phenotype.assertions.tsv.gz"
    # This is a builder output during fresh installation. Runtime consumes the
    # explicit key, with the legacy source-relative path only when it is unset.
    selected_hpo = config.get("hpo_assertions")
    checked_path(str(expected) if selected_hpo in (None, "", "null") else selected_hpo, "hpo_assertions")
    baseline = clinvar_stats_path(config, source_root)
    bundled = source_root / "data/ClinVar/clinvar_2star_stats.pkl.gz"
    # A fresh installer can recreate this cache from official ClinVar/MANE.
    # Defer only the genuine default: an explicitly selected external input or
    # a compatibility symlink must not be silently replaced with a new dataset.
    missing_default = (
        baseline == bundled
        and baseline.resolve() == bundled.absolute()
        and (not baseline.exists() or (baseline.is_file() and baseline.stat().st_size == 0))
    )
    if not (allow_missing_bundled_clinvar and missing_default):
        validate_selected_clinvar_stats(config, source_root)


def validate_bundled_snapshots(config: dict, source_root: Path) -> None:
    """The checkout supplies fixed scientific tables: existence/nonzero only.

    No unpickling, gzip/TSV parsing, digest, timestamp or upstream refresh is
    appropriate for these bundled snapshots. Their ordinary scientific reads
    remain in the original pipeline, in its original order.
    """
    assembly = {"hg19": "hg19", "GRCh37": "hg19", "hg38": "hg38", "GRCh38": "hg38"}.get(config.get("assembly"))
    if assembly is None:
        raise ValueError("Select hg19/GRCh37 or hg38/GRCh38")
    repeat_name = config.get("repeat_region_file_name") or f"ucsc.{assembly}.fasta.repeat.bed.gz"
    if (not isinstance(repeat_name, str) or repeat_name in (".", "..")
            or not re.fullmatch(r"[^/\\\r\n]+", repeat_name)):
        raise ValueError("repeat_region_file_name must be one filename beneath data/repeats, not an absolute or relative path")
    paths = [
        source_root / f"data/HCSeeker/HC_spots.{assembly}.tsv",
        source_root / f"data/gnomAD_RMC/gnomad_v2.1.1_rmc_flat.{assembly}.tsv",
        source_root / "data/alphamissense/alphamissense_mean_score.tsv",
        source_root / "data/vep/csq_interpretations.tsv",
        source_root / "data/repeats" / repeat_name,
    ]
    for key, relative in (
        ("hgnc_table", "data/hgnc/non_alt_loci_set.tsv"),
        ("dispensable_gene_list", "data/dispensable_genes/dispensable_gene_list.txt"),
        ("gofcards_mechanism_review_tsv", "data/gofcards/gofcards_mechanism_reviews.tsv"),
        ("mondo_disease_scope_overrides", "data/mondo/disease_scope_overrides.tsv"),
        ("gene_nonlof_mechanism_schema", "data/patho_mechanism/gene_nonlof_mechanism_curated_assertions.schema.json"),
    ):
        selected = config.get(key)
        path = source_root / relative if selected in (None, "", "null") else source_root / str(selected)
        if (key == "gofcards_mechanism_review_tsv" and selected not in (None, "", "null")
                and Path(os.path.abspath(path)) != Path(os.path.abspath(source_root / relative))
                and not path.exists() and not path.is_symlink()):
            # This explicit review is a build-only input. GoFCards' real
            # preflight requires it before any acquisition when rebuilding;
            # a current normalized cache does not need it for runtime reuse.
            continue
        paths.append(path)
    for path in paths:
        readable_file(str(path), "bundled snapshot or explicitly selected replacement")


def validate_fasta(value: str, assembly: str) -> None:
    key = f"gofcards_{assembly}_fasta"
    fasta = readable_file(value, key)
    index = readable_file(str(fasta) + ".fai", key + " index (run samtools faidx)")
    with index.open(encoding="utf-8") as handle:
        lengths = {row[0]: int(row[1]) for line in handle if len(row := line.rstrip().split("\t")) >= 2}
    if lengths.get("chr1", lengths.get("1")) != CHR1_LENGTH[assembly]:
        raise ValueError(f"{key} .fai does not identify the required {assembly} chromosome 1 length: {value}")


def validate_forward_chain(value: str) -> None:
    path = readable_file(value, "gofcards_hg19_to_hg38_chain")
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt", encoding="utf-8") as handle:
        # UCSC's target is the source coordinate system consumed by CrossMap.
        header = next((line.split() for line in handle if line.startswith("chain ")), [])
    if len(header) != 13 or header[2].removeprefix("chr") != "1" or header[7].removeprefix("chr") != "1" or (
            int(header[3]), int(header[8])) != (CHR1_LENGTH["hg19"], CHR1_LENGTH["hg38"]):
        raise ValueError(f"gofcards_hg19_to_hg38_chain must be the forward UCSC hg19-to-hg38 chain: {value}")


def validate_merged_cache(directory: str, version: str, assembly: str, *, allow_missing: bool = False) -> None:
    if not str(version).isdigit():
        raise ValueError("Set vep_installed_version to the selected numeric VEP release")
    root = checked_path(directory, "vep_cache_dir")
    selected = root / "homo_sapiens_merged" / f"{version}_{assembly}"
    cache = checked_path(str(selected), "GoFCards merged VEP cache")
    for parent in (Path(directory), root / "homo_sapiens_merged"):
        if os.path.lexists(parent) and not parent.is_dir():
            raise ValueError(f"GoFCards merged VEP cache parent is not a directory: {parent}")
    if allow_missing and not os.path.lexists(selected):
        return
    for relative in ("info.txt", "chr_synonyms.txt", "1/all_vars.gz"):
        readable_file(str(cache / relative), f"GoFCards {assembly} merged VEP {version} cache ({relative})")
    indexes = [checked_path(str(cache / ("1/all_vars.gz" + suffix)), "GoFCards merged VEP index")
               for suffix in (".csi", ".tbi")]
    if not any(path.is_file() and os.access(path, os.R_OK) and path.stat().st_size for path in indexes):
        raise ValueError(f"GoFCards {assembly} merged VEP {version} cache lacks a converted variant index: {cache}")
    metadata = {parts[0]: parts[1].strip() for line in (cache / "info.txt").read_text().splitlines()
                if len(parts := line.split(None, 1)) == 2}
    if any(metadata.get(key) != expected for key, expected in
           (("assembly", assembly), ("species", "homo_sapiens"), ("var_type", "tabix"))):
        raise ValueError(f"GoFCards {assembly} merged VEP {version} cache has mismatched assembly/species/converted metadata: {cache}")


def validate_gofcards(config: dict, *, main_preflight: bool = False,
                      allow_pinned_reference_defaults: bool = False,
                      allow_missing_merged_caches: bool = False,
                      source_root: Path | None = None) -> None:
    selected = {"hg19": "hg19", "GRCh37": "hg19", "hg38": "hg38", "GRCh38": "hg38"}.get(config.get("assembly"))
    if allow_pinned_reference_defaults:
        root = source_root or Path(__file__).resolve().parents[1]
        resources, _ = REFERENCES.resource_plan(config, root / "data", "gofcards", use_cache_root=True)
        # This resolves/validates every explicit reference without creating
        # directories, indexes or downloads. Only pinned missing defaults wait
        # for preparation; an original DRAGEN path is never a UCSC destination.
        REFERENCES.preflight(resources)
    else:
        for assembly in ("hg19", "hg38"):
            value = config.get(f"gofcards_{assembly}_fasta")
            if main_preflight and assembly == selected and value == config.get("ref_genome"):
                fasta = checked_path(value, f"gofcards_{assembly}_fasta")
                index = checked_path(str(fasta) + ".fai", f"gofcards_{assembly}_fasta index")
                if fasta.is_file() and index.is_file():
                    validate_fasta(value, assembly)
            else:
                validate_fasta(value, assembly)
        if config.get("gofcards_hg19_fasta") == config.get("gofcards_hg38_fasta"):
            raise ValueError("GoFCards hg19 and hg38 FASTAs must be distinct assembly-correct inputs")
        validate_forward_chain(config.get("gofcards_hg19_to_hg38_chain"))
    directory = config.get("vep_cache_dir")
    if main_preflight and config.get("cache_root"):
        folders = config.get("cache_folders") or {}
        directory = str(Path(config["cache_root"]) / folders.get("vep_cache", "VEP_cache"))
    for vep_assembly in ("GRCh37", "GRCh38"):
        validate_merged_cache(directory, str(config.get("vep_installed_version") or ""), vep_assembly,
                              allow_missing=allow_missing_merged_caches and os.environ.get("PRIVA_SMOKE_SKIP_VEP") != "1")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("config", type=Path)
    parser.add_argument("--source-root", type=Path, default=Path(__file__).resolve().parents[1])
    parser.add_argument("--stage", choices=("main", "manual-inputs", "snapshots", "mapping", "mapping-target", "mapping-source", "clinvar", "gofcards", "gofcards-main"), default="main")
    parser.add_argument("--prompt", action="store_true",
                        help="Manual-inputs stage only: interactively collect unset paths and save them together")
    parser.add_argument("--print-mapping", action="store_true", help="Print the validated resolved mapping path")
    parser.add_argument("--allow-missing-bundled-mapping", action="store_true",
                        help="Main preflight only: defer the one assembly-matched bundled mapping to its pinned installer")
    parser.add_argument("--allow-missing-bundled-clinvar", action="store_true",
                        help="Main preflight only: recreate an absent/empty default ClinVar gene-stat cache")
    parser.add_argument("--allow-pinned-reference-defaults", action="store_true",
                        help="GoFCards preflight only: allow its pinned missing reference/index/chain defaults to be installed")
    parser.add_argument("--allow-missing-merged-caches", action="store_true",
                        help="GoFCards installer preflight only: defer absent official caches to its cache preparation")
    parser.add_argument("--vep-version", help="Selected installed release, without changing configuration")
    args = parser.parse_args()
    try:
        if args.prompt and args.stage != "manual-inputs":
            raise ValueError("--prompt is restricted to --stage manual-inputs")
        if args.allow_missing_bundled_clinvar and args.stage != "main":
            raise ValueError("--allow-missing-bundled-clinvar is restricted to main preflight")
        if args.allow_pinned_reference_defaults and args.stage not in ("gofcards", "gofcards-main"):
            raise ValueError("--allow-pinned-reference-defaults is restricted to GoFCards stages")
        if args.allow_missing_merged_caches and args.stage not in ("gofcards", "gofcards-main"):
            raise ValueError("--allow-missing-merged-caches is restricted to GoFCards stages")
        if args.stage == "manual-inputs":
            manual_inputs(args.config, prompt=args.prompt)
            return 0
        with args.config.open(encoding="utf-8") as handle:
            config = yaml.safe_load(handle)
        if not isinstance(config, dict):
            raise ValueError("The configuration must be a YAML mapping")
        if args.vep_version is not None:
            config["vep_installed_version"] = args.vep_version
        if args.stage == "main":
            validate_main(config, args.source_root,
                          allow_missing_bundled_mapping=args.allow_missing_bundled_mapping,
                          allow_missing_bundled_clinvar=args.allow_missing_bundled_clinvar)
        elif args.stage == "snapshots":
            validate_bundled_snapshots(config, args.source_root)
        elif args.stage == "mapping":
            path = validate_mapping(config, args.source_root)
            if args.print_mapping:
                print(path)
        elif args.stage == "mapping-target":
            print(nm_to_enst_path(config, args.source_root))
        elif args.stage == "mapping-source":
            print(mapping_fetch_source(config, args.source_root))
        elif args.stage == "clinvar":
            validate_selected_clinvar_stats(config, args.source_root)
        else:
            validate_gofcards(config, main_preflight=args.stage == "gofcards-main",
                              allow_pinned_reference_defaults=args.allow_pinned_reference_defaults,
                              allow_missing_merged_caches=args.allow_missing_merged_caches,
                              source_root=args.source_root)
    except (OSError, ValueError, TypeError, yaml.YAMLError, BIOMART.ET.ParseError) as error:
        print(f"ERROR: {error}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
