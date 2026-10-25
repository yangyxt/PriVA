#!/usr/bin/env python3
"""Refresh ClinVar gene statistics after three calendar months of source-data age.

Entry points (run in the prepared PriVA environment):
  python scripts/clinvar_gene_refresh.py plan --config CONFIG --source-root ROOT
  python scripts/clinvar_gene_refresh.py install --config CONFIG --source-root ROOT
  python scripts/clinvar_gene_refresh.py acquire --plan PLAN --output INPUT_RECEIPT
  python scripts/clinvar_gene_refresh.py file-date --input ANNOTATED_VCF
  python scripts/clinvar_gene_refresh.py candidate-route --source-root ROOT
  python scripts/clinvar_gene_refresh.py candidate-begin --config CONFIG --source-root ROOT --input-vcf VCF
  python scripts/clinvar_gene_refresh.py candidate-complete --config CONFIG --source-root ROOT --input-vcf VCF

The baseline is checked only for existence/nonzero size and its UTC filesystem
mtime. Scientific parsing remains in the original installer/generator. Downloads,
working configurations, new statistics and receipts live in managed cache storage.
The original configuration is selected only after all ClinVar stages succeed.
"""
from __future__ import annotations

import argparse
import calendar
from contextlib import contextmanager
from datetime import date, datetime, timezone
import fcntl
import gzip
import hashlib
import json
import os
from pathlib import Path
import re
import shutil
import stat
import subprocess
import sys
import tempfile
import urllib.request
import zlib

import yaml

from installer_config_paths import UniqueLoader, render_config
from installer_prerequisites import clinvar_stats_path, validate_main


NCBI = "https://ftp.ncbi.nlm.nih.gov"
SUMMARY_URL = NCBI + "/pub/clinvar/tab_delimited/variant_summary.txt.gz"
MANE_URL = NCBI + "/refseq/MANE/MANE_human/release_1.4/MANE.GRCh38.v1.4.summary.txt.gz"
RECEIPT_KIND = "priva_clinvar_gene_statistics_refresh_v1"
CANDIDATE_KIND = "priva_clinvar_annotation_candidate_v1"
TIMEOUT = 30
HEADER_BYTES = 2 * 1024 * 1024
COMPRESSED_HEADER_BYTES = 1024 * 1024
EMPTY = (None, "", "null")
PROFILE_KEYS = (
    "assembly", "vep_installed_version", "ref_genome", "vep_cache_dir", "vep_plugins_dir",
    "vep_plugins_cachedir", "utr_annotator_file", "loeuf_prescore", "alphamissense_prescore",
    "spliceai_snv_prescore", "spliceai_indel_prescore", "primateai_prescore", "splicevault_prescore",
    "conservation_file", "loftee_repo", "human_ancestor_fasta", "loftee_conservation_file",
    "gerp_bigwig", "mavedb_file", "gnomad_vcf_chrX", "pext_enabled", "pext_mean_bw",
    "pext_tissue_bw_dir", "pext_tissues",
)
SCIENTIFIC_SOURCES = (
    "clinvar_stat_variants.py", "combine_annotations.py", "annotation_vcf.sh",
    "clinvar_annotation_contract.py", "common_bash_utils.sh", "stat_aachange_clinvar.py", "install_utils.sh", "hpo_penetrance.py",
)


def utc_now() -> str:
    return datetime.now(timezone.utc).isoformat()


def add_calendar_months(day: date, months: int = 3) -> date:
    """Clamp month ends; compare UTC calendar dates because fileDate has no time."""
    ordinal = day.year * 12 + day.month - 1 + months
    year, month = divmod(ordinal, 12)
    month += 1
    return date(year, month, min(day.day, calendar.monthrange(year, month)[1]))


def nonempty(path: Path, label: str) -> os.stat_result:
    info = path.stat()
    if not info.st_size:
        raise ValueError(f"Missing or empty {label}: {path}")
    return info


def digest(path: Path, algorithm: str = "sha256") -> str:
    value = hashlib.new(algorithm)
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            value.update(block)
    return value.hexdigest()


def identity(path: Path, **extra) -> dict:
    info = nonempty(path, "completed refresh artifact")
    return dict(extra, path=str(path), size=info.st_size, sha256=digest(path))


def baseline_metadata(path: Path) -> dict:
    """Inspect the bundled baseline without reading or interpreting its bytes."""
    value = {"path": str(path), "exists": path.exists()}
    if not path.exists():
        return value
    if not path.is_file():
        raise ValueError(f"Bundled ClinVar gene statistics is not a file: {path}")
    info = path.stat()
    value.update(size=info.st_size, mtime_ns=info.st_mtime_ns,
                 mtime_utc=datetime.fromtimestamp(info.st_mtime, timezone.utc).isoformat())
    return value


def atomic_text(path: Path, value: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    descriptor, temporary = tempfile.mkstemp(prefix=path.name + ".partial.", dir=path.parent)
    try:
        with os.fdopen(descriptor, "w") as handle:
            handle.write(value)
        if path.exists():
            os.chmod(temporary, path.stat().st_mode & 0o777)
        os.replace(temporary, path)
    finally:
        if os.path.exists(temporary):
            os.unlink(temporary)


def write_json(path: Path, value: dict) -> None:
    atomic_text(path, json.dumps(value, indent=2, sort_keys=True) + "\n")


def load_config(path: Path) -> dict:
    config = yaml.load(path.read_text(), Loader=UniqueLoader)
    if not isinstance(config, dict):
        raise ValueError("The configuration must be a YAML mapping")
    return config


def resource_policy(source_root: Path) -> dict:
    """Read active source/refresh choices from the shared resource registry."""
    path = source_root / "cache_resources.json"
    try:
        value = json.loads(path.read_text())
        if not isinstance(value, dict) or value.get("schema_version") != 1:
            raise ValueError("expected schema_version 1")
        resources = value["resources"]
        refresh_policy = resources["clinvar_gene_statistics"]["refresh"]
        months = refresh_policy["threshold_months"]
        automatic = refresh_policy["automatic"]
        if not isinstance(automatic, bool) or refresh_policy.get("policy") != "latest_filedate_vs_bundled_mtime":
            raise ValueError("expected boolean automatic and policy latest_filedate_vs_bundled_mtime")
        latest = resources["clinvar_vcf"]["source"]["latest_url_template"]
        summary = resources["clinvar_variant_summary"]["source"]["url"]
        mane = resources["mane_select"]["source"]["url"]
        if isinstance(months, bool) or not isinstance(months, int) or months <= 0:
            raise ValueError("threshold_months must be a positive integer")
        if (not isinstance(latest, str) or "{assembly}" not in latest
                or re.findall(r"\{([^}]+)\}", latest) != ["assembly"]
                or not latest.endswith("/clinvar.vcf.gz")):
            raise ValueError("latest_url_template must contain {assembly} and end in /clinvar.vcf.gz")
        for url in (latest, summary, mane):
            if not isinstance(url, str) or not url.startswith(NCBI + "/"):
                raise ValueError("ClinVar/MANE source URLs must use official NCBI HTTPS")
    except (OSError, ValueError, KeyError, TypeError) as error:
        raise ValueError(f"Cannot load ClinVar refresh policy from {path}: {error}") from error
    return {"registry": str(path), "registry_sha256": digest(path), "threshold_months": months, "automatic": automatic,
            "latest_url_template": latest, "summary_url": summary, "mane_url": mane}


def open_url(url: str, *, header_only: bool = False):
    headers = {"User-Agent": "PriVA-ClinVar-refresh/1", "Accept-Encoding": "identity"}
    if header_only:
        headers["Range"] = f"bytes=0-{COMPRESSED_HEADER_BYTES - 1}"
    return urllib.request.urlopen(urllib.request.Request(url, headers=headers), timeout=TIMEOUT)


class BoundedReader:
    """Bound compressed network reads even when an HTTP server ignores Range."""
    def __init__(self, source):
        self.source = source
        self.total = 0

    def read(self, count: int = -1) -> bytes:
        remaining = COMPRESSED_HEADER_BYTES - self.total
        count = min(count if count >= 0 else remaining + 1, remaining + 1)
        value = self.source.read(count)
        self.total += len(value)
        if self.total > COMPRESSED_HEADER_BYTES:
            raise ValueError("ClinVar compressed header exceeds the bounded metadata limit")
        return value


def header_date(handle) -> date:
    """Read the complete bounded header, without scanning variants or file mtimes."""
    result = None
    consumed = 0
    while consumed < HEADER_BYTES:
        line = handle.readline(min(65537, HEADER_BYTES - consumed + 1))
        consumed += len(line)
        if len(line) > 65536 or consumed > HEADER_BYTES:
            raise ValueError("ClinVar header exceeds the bounded metadata limit")
        if not line:
            break
        if line.startswith(b"##fileDate="):
            match = re.fullmatch(rb"##fileDate=(\d{8}|\d{4}-\d{2}-\d{2})\r?\n?", line)
            if not match:
                raise ValueError("ClinVar fileDate must be a YYYYMMDD or YYYY-MM-DD data date")
            stamp = match[1].decode()
            found = datetime.strptime(stamp, "%Y-%m-%d" if "-" in stamp else "%Y%m%d").date()
            if result is not None and result != found:
                raise ValueError("Conflicting ClinVar fileDate header values")
            result = found
        elif line.startswith(b"#CHROM\t"):
            if result is None:
                raise ValueError("ClinVar VCF header has no fileDate data date")
            return result
        elif not line.startswith(b"##"):
            raise ValueError("Expected a ClinVar VCF header before variant records")
    raise ValueError("ClinVar metadata ended before a complete VCF header/fileDate")


def local_file_date(path: Path) -> date:
    nonempty(path, "selected ClinVar VCF")
    with path.open("rb") as raw:
        magic = raw.read(2)
    opener = gzip.open if magic == b"\x1f\x8b" else open
    with opener(path, "rb") as handle:
        return header_date(handle)


def latest_metadata(assembly: str, *, url_template: str) -> dict:
    url = url_template.format(assembly=assembly)
    try:
        before_md5 = publisher_md5(url)
        with open_url(url, header_only=True) as response:
            with gzip.GzipFile(fileobj=BoundedReader(response), mode="rb") as handle:
                found = header_date(handle)
        if publisher_md5(url) != before_md5:
            raise ValueError("publisher MD5 changed during bounded header detection")
    except (OSError, EOFError, ValueError, zlib.error) as error:
        raise ValueError(f"Cannot determine latest official ClinVar fileDate from {url}: {error}") from error
    return {"url": url, "file_date": found.isoformat(), "observed_utc": utc_now(),
            "publisher_md5": before_md5,
            "date_source": "VCF ##fileDate; not download or local cache mtime"}


def receipt_path(output: Path) -> Path:
    return Path(str(output) + ".ready.json")


def owned_receipt(output: Path, *, assembly: str, vep_major: str) -> dict | None:
    path = receipt_path(output)
    if not path.exists():
        return None
    value = json.loads(path.read_text())
    if value.get("kind") != RECEIPT_KIND:
        return None
    if (value.get("status") != "complete" or value.get("assembly") != assembly
            or value.get("vep_major") != vep_major
            or Path(value["output"]["path"]).absolute() != output.absolute()):
        raise ValueError(f"Completed ClinVar refresh receipt does not match the selected profile: {path}")
    current = identity(output)
    if current != value["output"]:
        raise ValueError(f"Generated ClinVar statistics differ from their completed receipt: {output}")
    return value


def selected_vcf(config: dict, root: Path, ucsc: str) -> Path:
    value = config.get("clinvar_vcf")
    if value in EMPTY:
        return root / f"clinvar.{ucsc}.vep.vcf.gz"
    path = Path(value)
    return path if path.is_absolute() or "/" in str(value) else root / path


def scientific_profile(config: dict, source_root: Path, policy: dict) -> dict:
    """Bind scientific paths/code and actual input metadata, without hashing bundles.

    Large scientific resources use path/size/mtime plus their small index/cache
    metadata. Code, plugins and small external inputs also have SHA256 identities.
    Bundled data are never opened here for content validation.
    """
    values = {key: config.get(key) for key in PROFILE_KEYS}
    if str(config.get("pext_enabled", "false")).lower() != "true":
        for key in ("pext_mean_bw", "pext_tissue_bw_dir", "pext_tissues"):
            values.pop(key, None)
    files = {}
    def record(path: Path, *, code: bool = False) -> None:
        path = path.absolute()
        resolved = path.resolve()
        value = {"path": str(path), "resolved": str(resolved), "exists": path.exists()}
        if path.exists():
            info = path.stat()
            value["directory"] = path.is_dir()
            # Shared cache/plugin containers gain unrelated files during a
            # normal install. Bind actual selected resources, not directory mtimes.
            if path.is_file():
                value.update(size=info.st_size, mtime_ns=info.st_mtime_ns)
            bundled = resolved.is_relative_to(source_root.resolve() / "data")
            if path.is_file() and (code or (not bundled and info.st_size <= 1024 * 1024)):
                value["sha256"] = digest(path)
        files[str(path)] = value
    for key, value in values.items():
        if key in ("assembly", "vep_installed_version", "pext_enabled", "pext_tissues") or value in EMPTY:
            continue
        if isinstance(value, str):
            path = Path(value)
            if key == "loeuf_prescore":
                # Match the existing VEP shell's bundled basename resolution.
                path = source_root / "data/loeuf" / path.name
            record(path)
            if path.is_file():
                for suffix in (".tbi", ".csi", ".fai"):
                    index = Path(str(path) + suffix)
                    if index.exists():
                        record(index)
    for name in SCIENTIFIC_SOURCES:
        record(source_root / "scripts" / name, code=True)
    plugins = config.get("vep_plugins_dir")
    if plugins not in EMPTY:
        for name in ("UTRAnnotator", "LOEUF", "AlphaMissense", "SpliceAI", "PrimateAI", "SpliceVault", "Conservation", "NMD", "MaveDB"):
            record(Path(plugins) / (name + ".pm"), code=True)
    if config.get("loftee_repo") not in EMPTY:
        for path in sorted(Path(config["loftee_repo"]).glob("*.p[ml]")):
            record(path, code=True)
    if config.get("gnomad_vcf_chrX") not in EMPTY:
        for chrom in [str(value) for value in range(1, 23)] + ["X", "Y"]:
            path = Path(str(config["gnomad_vcf_chrX"]).replace(".chrX.vcf", ".chr" + chrom + ".vcf"))
            record(path)
            record(Path(str(path) + ".tbi"))
    if config.get("vep_cache_dir") not in EMPTY:
        assembly = "GRCh37" if config["assembly"] in ("hg19", "GRCh37") else "GRCh38"
        folder = Path(config["vep_cache_dir"]) / "homo_sapiens_merged" / f"{config['vep_installed_version']}_{assembly}"
        for name in ("info.txt", "chr_synonyms.txt"):
            record(folder / name)
    if str(config.get("pext_enabled", "false")).lower() == "true" and config.get("pext_tissue_bw_dir") not in EMPTY:
        tissues = config.get("pext_tissues") or []
        if isinstance(tissues, str):
            tissues = yaml.safe_load(tissues)
            if isinstance(tissues, str):
                tissues = [name.strip() for name in tissues.split(",") if name.strip()]
        if isinstance(tissues, list):
            for tissue in tissues:
                record(Path(config["pext_tissue_bw_dir"]) / "tissues" / (str(tissue) + ".bw"))
    executable = shutil.which("vep")
    if executable:
        record(Path(executable), code=True)
    package_names = ("python-", "pandas-", "numpy-", "pysam-", "pyarrow-", "ensembl-vep-", "htslib-", "bcftools-", "samtools-", "perl-")
    for path in sorted((Path(sys.prefix) / "conda-meta").glob("*.json")):
        if path.name.startswith(package_names):
            record(path, code=True)
    return {"config": values, "files": files,
            "python": {"executable": sys.executable, "version": sys.version},
            "source_policy": {key: policy[key] for key in ("latest_url_template", "summary_url", "mane_url", "automatic", "threshold_months")}}


def profile_signature(profile: dict | None) -> str:
    """Name an existing profile deterministically, without observation timestamps.

    A recorded content hash supersedes mtime for code/small inputs. Large inputs
    retain the metadata identity already used by the existing refresh receipt;
    this procedure never reads those payloads or hashes bundled data.
    """
    if not isinstance(profile, dict):
        return ""
    stable = json.loads(json.dumps(profile))
    for value in stable.get("files", {}).values():
        if "sha256" in value:
            value.pop("mtime_ns", None)
    return hashlib.sha256(json.dumps(stable, sort_keys=True, separators=(",", ":")).encode()).hexdigest()


def candidate_route(source_root: Path) -> dict:
    """Check source wiring before cache paths exist; never qualify a profile."""
    script = source_root / "scripts/annotation_vcf.sh"
    text = script.read_text()
    functions = ("anno_VEP_data", "anno_VEP_cache_candidate", "_priva_anno_VEP_data",
                 "priva_require_legacy_spliceai")
    for name in functions:
        if not re.search(r"^function " + re.escape(name) + r"\(\) \{", text, re.M):
            raise ValueError(f"ClinVar candidate preparation entry point is missing: {name}")
    guard = re.search(r"^function priva_require_legacy_spliceai\(\) \{.*?^\}", text, re.M | re.S)
    return {"kind": CANDIDATE_KIND, "status": "preparation_route_available_inputs_not_checked",
            "scientific_qualified": False, "source": identity(script),
            "helper_source": identity(source_root / "scripts/clinvar_gene_refresh.py"),
            "runtime_guard_sha256": hashlib.sha256(guard[0].encode()).hexdigest()}


def candidate_receipt_path(vcf: Path, *, pending: bool = False) -> Path:
    return Path(str(vcf) + (".installation-candidate.pending.json" if pending
                            else ".installation-candidate.ready.json"))


def candidate_context(config_path: Path, source_root: Path, input_vcf: Path) -> dict:
    """Bind a managed preparation to its actual inputs, without scientific activation.

    Large SpliceAI2 bodies have already been verified by their adapter. Reuse its
    bounded publication checks and receipt instead of hashing those caches again.
    """
    from spliceai_cache_input import _completed_projection, validate_pair

    config = load_config(config_path)
    root_value = config.get("cache_root")
    if root_value in EMPTY or not Path(root_value).is_absolute():
        raise ValueError("ClinVar candidate preparation requires an absolute configured cache_root")
    root = Path(root_value).resolve()
    assembly = {"hg19": "GRCh37", "GRCh37": "GRCh37", "hg38": "GRCh38", "GRCh38": "GRCh38"}.get(config.get("assembly"))
    if assembly is None:
        raise ValueError("ClinVar candidate preparation requires a configured assembly")
    ucsc = "hg19" if assembly == "GRCh37" else "hg38"
    folder = Path(config["clinvar_vcf_dir"])
    selected = selected_vcf(config, folder, ucsc)
    vcf = input_vcf.absolute()
    if (not folder.is_absolute() or not folder.resolve().is_relative_to(root)
            or not vcf.is_relative_to(Path(root_value).absolute())
            or not vcf.resolve().is_relative_to(folder.resolve())
            or vcf.resolve() != selected.resolve()):
        raise ValueError("ClinVar candidate output must be the selected managed VCF inside cache_root/clinvar_vcf_dir")
    # The original shell replaces the first '.vcf' suffix; require its output to
    # be exactly this owned file rather than allow an untracked alternate path.
    if not str(vcf).endswith(".vcf.gz") or str(vcf).split(".vcf", 1)[0] + ".vcf.gz" != str(vcf):
        raise ValueError("ClinVar candidate preparation requires a managed .vcf.gz path without earlier .vcf components")
    nonempty(vcf, "managed ClinVar candidate input")
    snv, indel = (config.get(key) for key in ("spliceai_snv_prescore", "spliceai_indel_prescore"))
    if snv in EMPTY or indel in EMPTY:
        raise ValueError("ClinVar candidate preparation requires both selected SpliceAI VCFs")
    pair = validate_pair(assembly, snv, indel)
    if pair["format"] not in ("vcf_v1", "vcf_spliceai2"):
        raise ValueError("ClinVar candidate preparation requires installed SpliceAI VCFs, not native publisher TSVs")
    adapter = None
    if pair["format"] == "vcf_spliceai2":
        with _completed_projection(snv, indel):
            adapter = identity(Path(snv).resolve().parent / "adapter.ready.json")
    return {"kind": CANDIDATE_KIND, "config_path": str(config_path.absolute()),
            "config_sha256": digest(config_path), "source_root": str(source_root.absolute()),
            "input_vcf": str(vcf), "assembly": assembly, "spliceai_format": pair["format"],
            "adapter_receipt": adapter, "route": candidate_route(source_root),
            "scientific_profile": scientific_profile(config, source_root, resource_policy(source_root)),
            "scientific_qualified": False}


def candidate_begin(config_path: Path, source_root: Path, input_vcf: Path) -> dict:
    value = candidate_context(config_path, source_root, input_vcf)
    vcf = Path(value["input_vcf"])
    info = vcf.stat()
    value.update(status="preparing_annotation_only", started_utc=utc_now(),
                 input_before={"path": str(vcf), "size": info.st_size, "mtime_ns": info.st_mtime_ns})
    write_json(candidate_receipt_path(vcf, pending=True), value)
    return value


def candidate_complete(config_path: Path, source_root: Path, input_vcf: Path) -> dict:
    from clinvar_annotation_contract import validate_header

    vcf = input_vcf.absolute()
    pending = candidate_receipt_path(vcf, pending=True)
    before = json.loads(pending.read_text())
    current = candidate_context(config_path, source_root, vcf)
    if (before.get("status") != "preparing_annotation_only"
            or any(before.get(key) != value for key, value in current.items())):
        raise ValueError("ClinVar candidate configuration/scientific inputs changed during preparation")
    contract = validate_header(vcf, "vep")
    value = dict(before, status="prepared_annotation_only", completed_utc=utc_now(),
                 stage="vep_before_gnomad_enrichment",
                 output=identity(vcf), output_index=identity(Path(str(vcf) + ".tbi")),
                 annotation_contract=contract,
                 qualification="installation candidate only; ACMG/ClinVar/CADD and public/family qualification remain pending")
    write_json(candidate_receipt_path(vcf), value)
    pending.unlink()
    return value


def validate_reuse(prior: dict, plan: dict) -> None:
    """A ready statistic alone cannot authorize stale annotation/input reuse."""
    if profile_signature(prior.get("scientific_profile")) != profile_signature(plan["scientific_profile"]):
        raise ValueError("ClinVar refresh scientific configuration/input/source identity changed; "
                         "the previous statistic was not selected for this different profile")
    if prior["latest"].get("publisher_md5") != plan["latest"].get("publisher_md5"):
        raise ValueError("Official ClinVar bytes changed within the recorded fileDate; no stale result was selected")
    for label, recorded in prior.get("completed_artifacts", {}).items():
        try:
            current = identity(Path(recorded["path"]))
        except OSError as error:
            raise ValueError(f"Completed ClinVar artifact is missing: {label}: {recorded['path']}") from error
        if current != recorded:
            raise ValueError(f"Completed ClinVar artifact changed: {label}: {recorded['path']}")
    if set(prior.get("completed_artifacts", {})) != {"vcf", "vcf_index", "aa", "splice"}:
        raise ValueError("ClinVar ready receipt does not bind all completed annotation/AA/splice artifacts")
    for recorded in prior["inputs"]["sources"].values():
        current = identity(Path(recorded["path"]))
        if any(current[key] != recorded[key] for key in ("path", "size", "sha256")):
            raise ValueError("A pinned ClinVar/MANE regeneration input changed; no stale result was selected")
    if plan["explicit_vcf"] and identity(Path(plan["input_vcf"])) != prior["completed_artifacts"]["vcf"]:
        raise ValueError("Explicit ClinVar VCF identity differs from the completed refresh, even though fileDate matches; "
                         "the explicit selection was not replaced")


def plan_refresh(config_path: Path, source_root: Path) -> dict:
    """Inspect only bundled metadata; probe NCBI for a default/owned selection."""
    config = load_config(config_path)
    pair = {"hg19": ("GRCh37", "hg19"), "GRCh37": ("GRCh37", "hg19"),
            "hg38": ("GRCh38", "hg38"), "GRCh38": ("GRCh38", "hg38")}.get(config.get("assembly"))
    if pair is None:
        raise ValueError("ClinVar refresh requires hg19/GRCh37 or hg38/GRCh38")
    assembly, ucsc = pair
    major = str(config.get("vep_installed_version", ""))
    if not re.fullmatch(r"[1-9][0-9]*", major):
        raise ValueError("ClinVar refresh requires the configured numeric VEP major")
    baseline = source_root / "data/ClinVar/clinvar_2star_stats.pkl.gz"
    selected = clinvar_stats_path(config, source_root)
    is_baseline = selected.absolute() == baseline.absolute()
    if not is_baseline:
        nonempty(selected, "selected ClinVar gene statistics")
    known = owned_receipt(selected, assembly=assembly, vep_major=major) if not is_baseline else None
    plan = {"kind": RECEIPT_KIND, "config_path": str(config_path), "source_root": str(source_root),
            "config_sha256": digest(config_path), "assembly": assembly, "ucsc_assembly": ucsc,
            "vep_major": major, "selected_stat": str(selected),
            "action": "preserve", "reason": "explicitly selected statistics are retained"}
    if not is_baseline and known is None:
        return plan
    policy = resource_policy(source_root)
    baseline_state = baseline_metadata(baseline)
    missing_baseline = not baseline_state.get("size")
    if not policy["automatic"] and not missing_baseline:
        plan.update(resource_policy=policy, reason="automatic ClinVar gene-statistics refresh is disabled in the resource registry")
        return plan
    boundary = None if missing_baseline else add_calendar_months(
        datetime.fromisoformat(baseline_state["mtime_utc"]).date(), policy["threshold_months"])
    plan.update(resource_policy=policy, baseline=baseline_state,
                baseline_missing=missing_baseline)
    if boundary is not None:
        plan["baseline"]["refresh_boundary_utc_date"] = boundary.isoformat()
    latest = latest_metadata(assembly, url_template=policy["latest_url_template"])
    plan["latest"] = latest
    if boundary is not None and date.fromisoformat(latest["file_date"]) < boundary:
        plan.update(reason=f"latest ClinVar data date is before the {policy['threshold_months']}-calendar-month boundary")
        return plan
    root = Path(config["clinvar_vcf_dir"])
    profile = scientific_profile(config, source_root, policy)
    signature = profile_signature(profile)
    release_dir = root / "gene_statistics" / assembly / latest["file_date"].replace("-", "")
    build = release_dir / ("vep" + major + "-" + signature[:16])
    output = build / "clinvar_2star_stats.pkl.gz"
    prior = owned_receipt(output, assembly=assembly, vep_major=major)
    # Reuse an existing unsuffixed completion only for its exact profile. A
    # different profile gets a separate workspace and cannot relabel old bytes.
    legacy_build = release_dir / ("vep" + major)
    legacy_output = legacy_build / "clinvar_2star_stats.pkl.gz"
    if prior is None and receipt_path(legacy_output).is_file():
        legacy_prior = owned_receipt(legacy_output, assembly=assembly, vep_major=major)
        if legacy_prior is not None and profile_signature(legacy_prior.get("scientific_profile")) == signature:
            build, output, prior = legacy_build, legacy_output, legacy_prior
    input_vcf = selected_vcf(config, root, ucsc)
    default_vcf = root / f"clinvar.{ucsc}.vep.vcf.gz"
    prior_vcf = (known or {}).get("config_updates", {}).get("clinvar_vcf")
    managed = input_vcf.absolute() == default_vcf.absolute() or (
        prior_vcf is not None and input_vcf.absolute() == Path(prior_vcf).absolute()
        and (known or {}).get("managed_root") == str(root))
    if not managed and local_file_date(input_vcf).isoformat() != latest["file_date"]:
        raise ValueError(f"Explicit ClinVar VCF {input_vcf} has fileDate {local_file_date(input_vcf)}; "
                         f"latest official ClinVar is {latest['file_date']}. No replacement was selected; "
                         "do not combine a newer summary with the explicitly selected older/different release.")
    plan.update(action="reuse" if prior is not None and prior["data_date"] == latest["file_date"] else "rebuild",
                reason=("missing/empty bundled ClinVar gene statistics require managed reconstruction" if missing_baseline
                        else f"latest ClinVar data date reaches the {policy['threshold_months']}-calendar-month boundary"),
                managed_root=str(root), build_dir=str(build), output_path=str(output),
                input_vcf=str(input_vcf), explicit_vcf=not managed)
    plan["explicit_vcf_identity"] = identity(input_vcf) if not managed else None
    plan["scientific_profile"] = profile
    plan["scientific_profile_sha256"] = signature
    if prior is not None and prior["data_date"] == latest["file_date"]:
        validate_reuse(prior, plan)
    return plan


def publisher_md5(url: str) -> str:
    with open_url(url + ".md5") as response:
        text = response.read(4097)
    if len(text) > 4096 or not re.match(rb"^[0-9a-fA-F]{32}(?:\s|$)", text):
        raise ValueError(f"Invalid publisher MD5 manifest: {url}.md5")
    return text.split()[0].decode().lower()


def acquire_file(url: str, destination: Path, *, expected_md5: str | None = None) -> dict:
    """Reuse only completed external bytes; publish downloads after their checks."""
    sidecar = Path(str(destination) + ".source.json")
    reusable = False
    recorded = json.loads(sidecar.read_text()) if sidecar.exists() else {}
    http_provenance = recorded.get("http_response") if recorded.get("url") == url else None
    if destination.exists():
        if expected_md5 is not None:
            reusable = digest(destination, "md5") == expected_md5
        elif sidecar.exists():
            reusable = recorded.get("url") == url and recorded.get("sha256") == digest(destination)
    if not reusable:
        destination.parent.mkdir(parents=True, exist_ok=True)
        descriptor, temporary = tempfile.mkstemp(prefix=destination.name + ".partial.", dir=destination.parent)
        try:
            with os.fdopen(descriptor, "wb") as output, open_url(url) as response:
                headers = getattr(response, "headers", {})
                http_provenance = {"last_modified": headers.get("Last-Modified"),
                                   "etag": headers.get("ETag"), "content_length": headers.get("Content-Length"),
                                   "retrieved_utc": utc_now()}
                shutil.copyfileobj(response, output, length=1024 * 1024)
            temporary_path = Path(temporary)
            nonempty(temporary_path, "downloaded ClinVar refresh input")
            if expected_md5 is not None and digest(temporary_path, "md5") != expected_md5:
                raise ValueError(f"Downloaded bytes differ from publisher MD5: {url}")
            if destination.name.endswith(".gz"):
                with gzip.open(temporary_path, "rb") as handle:
                    for _ in iter(lambda: handle.read(1024 * 1024), b""):
                        pass
            os.replace(temporary, destination)
        finally:
            if os.path.exists(temporary):
                os.unlink(temporary)
    result = identity(destination, url=url, md5=digest(destination, "md5"))
    if http_provenance is not None:
        result["http_response"] = http_provenance
    if expected_md5 is not None:
        result["publisher_md5"] = expected_md5
    write_json(sidecar, result)
    return result


def acquire_inputs(plan: dict) -> dict:
    folder = Path(plan["build_dir"]) / "inputs"
    folder.mkdir(parents=True, exist_ok=True)
    if plan["explicit_vcf"]:
        sources = {"vcf": identity(Path(plan["input_vcf"]), kind="explicit_complete_annotated_vcf")}
    else:
        stamp = plan["latest"]["file_date"].replace("-", "")
        # Weekly releases need not have a dated archive URL. Pin the moving
        # endpoint with publisher MD5 and verify it again after acquisition.
        url = plan["resource_policy"]["latest_url_template"].format(assembly=plan["assembly"])
        expected = publisher_md5(url)
        if expected != plan["latest"].get("publisher_md5"):
            raise ValueError("Official ClinVar changed after latest-data detection; no replacement selected")
        vcf = folder / f"clinvar_{stamp}.{plan['assembly']}.md5_{expected}.vcf.gz"
        sources = {"vcf": acquire_file(url, vcf, expected_md5=expected),
                   "index": acquire_file(url + ".tbi", Path(str(vcf) + ".tbi"))}
    if local_file_date(Path(sources["vcf"]["path"])).isoformat() != plan["latest"]["file_date"]:
        raise ValueError("Pinned ClinVar input fileDate differs from the detected latest release")
    summary_url = plan["resource_policy"]["summary_url"]
    mane_url = plan["resource_policy"]["mane_url"]
    summary_md5 = publisher_md5(summary_url)
    sources["summary"] = acquire_file(summary_url, folder / f"variant_summary.md5_{summary_md5}.txt.gz",
                                      expected_md5=summary_md5)
    sources["mane"] = acquire_file(mane_url, folder / mane_url.rsplit("/", 1)[1])
    if publisher_md5(summary_url) != summary_md5:
        raise ValueError("ClinVar variant_summary publisher MD5 changed during acquisition; no statistics selected")
    after = latest_metadata(plan["assembly"], url_template=plan["resource_policy"]["latest_url_template"])
    if (after["file_date"] != plan["latest"]["file_date"]
            or after["publisher_md5"] != plan["latest"]["publisher_md5"]):
        raise ValueError("Official ClinVar changed during acquisition; no partial release was selected")
    mane_release = re.search(r"/release_([^/]+)/", mane_url)
    value = {"kind": RECEIPT_KIND, "sources": sources, "latest": plan["latest"], "acquired_utc": utc_now(),
             "mane_release": mane_release[1] if mane_release else mane_url,
             "summary_release_identity": "independent publisher MD5 snapshot; not inferred from filename or VCF date"}
    write_json(Path(plan["build_dir"]) / "inputs.ready.json", value)
    return value


def run_original(source_root: Path, function: str, *arguments: str | Path) -> None:
    subprocess.run(["bash", str(source_root / "scripts/install_utils.sh"), function,
                    *(str(value) for value in arguments)], check=True)


@contextmanager
def build_lock(folder: Path):
    folder.mkdir(parents=True, exist_ok=True)
    with (folder / ".refresh.lock").open("a") as handle:
        fcntl.flock(handle, fcntl.LOCK_EX)
        yield




def deployment_context(config_path: Path) -> dict:
    """Resolve ownership without inferring a profile from old file contents."""
    config_path = config_path.absolute()
    config = load_config(config_path)
    assembly = {"hg19": "GRCh37", "GRCh37": "GRCh37",
                "hg38": "GRCh38", "GRCh38": "GRCh38"}.get(config.get("assembly"))
    major = str(config.get("vep_installed_version", ""))
    if assembly is None or not re.fullmatch(r"[1-9][0-9]*", major):
        raise ValueError("ClinVar deployment requires a configured assembly and numeric VEP major")
    folder = Path(config["clinvar_vcf_dir"])
    if not folder.is_absolute():
        raise ValueError("ClinVar deployment requires an absolute clinvar_vcf_dir")
    folder = Path(os.path.abspath(folder))
    ucsc = "hg19" if assembly == "GRCh37" else "hg38"
    target = Path(os.path.abspath(selected_vcf(config, folder, ucsc)))
    default = folder / f"clinvar.{ucsc}.vep.vcf.gz"
    owned = target == default
    if target.is_symlink():
        raise ValueError(f"Linked ClinVar annotation preserved without replacement: {target}")
    root_value = config.get("cache_root")
    root = Path(root_value).absolute() if isinstance(root_value, str) and root_value.startswith("/") else None
    if owned:
        if root is None or not folder.resolve().is_relative_to(root.resolve()):
            raise ValueError("Owned ClinVar output must remain inside the configured cache_root; prior paths preserved")
        if str(target).split(".vcf", 1)[0] + ".vcf.gz" != str(target):
            raise ValueError("Owned ClinVar output cannot have an earlier .vcf path component")
    return dict(config=config, config_path=config_path, folder=folder, target=target,
                root=root, assembly=assembly, ucsc=ucsc, major=major, owned=owned)


def deployment_header(context: dict, input_file: Path | None = None) -> dict:
    from alphamissense_cache import check_vep_annotation

    return check_vep_annotation(input_file=context["target"] if input_file is None else input_file, assembly=context["assembly"],
                                vep_major=context["major"], resource="ClinVar")


def retirement_members(context: dict) -> list[Path]:
    """Fixed names replaced by acquisition/VEP; old shards are handled later."""
    folder, target = context["folder"], context["target"]
    prefix = target.name.split(".vcf", 1)[0]
    payloads = [target,
                folder / f"clinvar.{context['assembly']}.vcf.gz",
                folder / f"clinvar.{context['ucsc']}.vcf.gz"]
    members = [Path(str(path) + suffix) for path in payloads for suffix in ("", ".tbi", ".csi")]
    members += [Path(str(path)[:-3]) for path in payloads]
    members += [Path(str(target) + suffix) for suffix in
                (".installation-candidate.pending.json", ".installation-candidate.ready.json",
                 ".gnomAD.publish.pending")]
    members += [folder / (prefix + suffix) for suffix in (".gnomAD.log", ".vep.stats.html")]
    return list(dict.fromkeys(members))


def shard_retirement_members(context: dict, chromosomes: list[str]) -> list[Path]:
    """Use the unchanged writer's literal selected names, never inferred contigs."""
    folder, target = context["folder"], context["target"]
    prefix = target.name.split(".vcf", 1)[0]
    members = []
    for chrom in chromosomes:
        if not chrom or any(character in chrom for character in "/\x00\r\n\t"):
            raise ValueError("Unsafe literal ClinVar shard name; all paths preserved")
        for ending in (".vcf.gz", ".gnomAD.vcf.gz"):
            shard = folder / f"{prefix}.{chrom}{ending}"
            members += [Path(str(shard) + suffix) for suffix in ("", ".tbi", ".csi")]
            members.append(Path(str(shard)[:-3]))
            if ending == ".gnomAD.vcf.gz":
                members.append(Path(str(shard) + ".gnomAD.publish.pending"))
    return list(dict.fromkeys(members))


def retirement_identity(path: Path) -> dict:
    info = path.lstat()
    if not stat.S_ISREG(info.st_mode):
        raise ValueError(f"Nonregular/linked ClinVar retirement member preserved: {path}")
    return {name: str(getattr(info, "st_" + name)) for name in
            ("dev", "ino", "mode", "nlink", "size", "mtime_ns", "ctime_ns")}


def persisted_identity_matches(recorded: dict, observed: dict) -> bool:
    """Compare a retry observation without the native node's device number."""
    # 2026-10-23 CT: the same shared file has different st_dev on native nodes.
    # Retain all seven raw fields; only a persisted comparison omits dev. Header
    # reads, inventory/move checks and checks within this invocation stay exact.
    fields = {"dev", "ino", "mode", "nlink", "size", "mtime_ns", "ctime_ns"}
    return (isinstance(recorded, dict) and isinstance(observed, dict)
            and set(recorded) == set(observed) == fields
            and all(isinstance(value[name], str) for value in (recorded, observed) for name in fields)
            and all(re.fullmatch(r"[0-9]+", value["dev"]) for value in (recorded, observed))
            and all(recorded[name] == observed[name] for name in fields - {"dev"}))


def record_retry_observation(pending: dict, *, protected_sources: dict,
                             rebuilt_target: dict | None = None) -> None:
    """Append current raw identities without rewriting historical before values."""
    observations = pending.setdefault("retry_observations", [])
    if not isinstance(observations, list):
        raise ValueError("Malformed ClinVar retry observations; previous identities preserved")
    value = dict(observed_utc=utc_now(), phase=pending["phase"], protected_sources=protected_sources)
    if rebuilt_target is not None:
        value["rebuilt_target"] = rebuilt_target
    observations.append(value)


def retirement_inventory(context: dict, source_vcf: Path | None = None, *, members=None,
                         protect_target: bool = False) -> list:
    """Reject foreign aliases before creating an archive or moving any member."""
    root = context["root"].resolve()
    protected = [context["config_path"]]
    if protect_target:
        protected += [context["target"], Path(str(context["target"]) + ".tbi"),
                      Path(str(context["target"]) + ".csi")]
    if source_vcf is not None:
        protected += [source_vcf.absolute(), Path(str(source_vcf.absolute()) + ".tbi"),
                      Path(str(source_vcf.absolute()) + ".csi")]
    moves = []
    for member in retirement_members(context) if members is None else members:
        if not member.parent.resolve().is_relative_to(root):
            raise ValueError(f"ClinVar retirement parent escapes cache_root; preserved: {member}")
        if not os.path.lexists(member):
            continue
        identity_before = retirement_identity(member)
        if any(member.resolve() == item.resolve() or (item.exists() and member.samefile(item))
               for item in protected):
            raise ValueError(f"ClinVar retirement would overwrite a selected config/source; preserved: {member}")
        moves.append((member, identity_before))
    return moves


def deployment_pending_path(context: dict) -> Path:
    return context["folder"].resolve() / ".clinvar-deploy-locks" / context["target"].name / "retirement.pending.json"


def deployment_pending(context: dict) -> dict | None:
    path = deployment_pending_path(context)
    if path.is_symlink():
        raise ValueError("Linked ClinVar recovery journal preserved")
    if not path.exists():
        return None
    value = json.loads(path.read_text())
    expected = dict(kind="priva_clinvar_annotation_recovery_v1", target=str(context["target"]),
                    config=str(context["config_path"]), assembly=context["assembly"], vep_major=context["major"])
    if (any(value.get(key) != item for key, item in expected.items())
            or value.get("phase") not in ("awaiting_vep", "retiring_shards", "shards_retained")):
        raise ValueError("ClinVar recovery journal belongs to a different target/configuration; preserved")
    return value


def pinned_source_identity(context: dict, source_vcf: Path | None) -> dict:
    """A supplied source cannot also occupy the controlled overwrite namespace."""
    result = {}
    if source_vcf is None:
        return result
    prefix = context["target"].name.split(".vcf", 1)[0] + "."
    for path in (source_vcf.absolute(), Path(str(source_vcf.absolute()) + ".tbi"),
                 Path(str(source_vcf.absolute()) + ".csi")):
        if not os.path.lexists(path):
            continue
        info = retirement_identity(path)
        if (int(info["nlink"]) != 1 or
                (path.resolve().parent == context["folder"].resolve() and
                 path.resolve().name.startswith(prefix))):
            raise ValueError("Pinned ClinVar source/index may alias a managed overwrite member; preserved: " + str(path))
        result[str(path)] = info
    return result


def retire_owned_cache(config_path: Path, *, reason: str, source_vcf: Path | None = None,
                       expected_target: dict | None = None, expected_members: dict | None = None) -> dict | None:
    """Retain exact overwrite members. Caller holds the distinct deployment lock.

    The journal is written before renames. Companions precede the main VCF, so
    interruption cannot make an older annotation eligible for a new warm reuse.
    Retry retains any remaining members in a new journal; no payload is deleted.
    """
    context = deployment_context(config_path)
    if not context["owned"]:
        raise ValueError("Explicit ClinVar annotations are preserved; automatic retirement is limited to the managed default")
    pending = deployment_pending(context)
    if pending is not None and pending["phase"] != "awaiting_vep":
        raise ValueError("ClinVar shard retirement is pending; the new VEP output was preserved")
    if os.path.lexists(context["target"]):
        from alphamissense_cache import AnnotationHeaderChangedError

        before = retirement_identity(context["target"])
        if expected_target is not None and before != expected_target:
            raise AnnotationHeaderChangedError("ClinVar annotation changed after its header check; preserved")
        try:
            deployment_header(context)
        except AnnotationHeaderChangedError:
            raise
        except (ValueError, EOFError, gzip.BadGzipFile, zlib.error):
            pass  # Only a stable, owned unreadable/mismatched header may retire.
        if retirement_identity(context["target"]) != before:
            raise AnnotationHeaderChangedError("ClinVar annotation changed after its header check; preserved")
    elif expected_target is not None:
        raise ValueError("ClinVar annotation disappeared after its header check; no retirement started")
    moves = retirement_inventory(context, source_vcf)
    protected_sources = pinned_source_identity(context, source_vcf)
    for path, expected in (expected_members or {}).items():
        if retirement_identity(path) != expected:
            raise ValueError("ClinVar repair source changed after its header check; no retirement started")
    if pending is None:
        pending = dict(kind="priva_clinvar_annotation_recovery_v1", phase="awaiting_vep",
                       target=str(context["target"]), config=str(context["config_path"]),
                       assembly=context["assembly"], vep_major=context["major"],
                       started_utc=utc_now(), scientific_qualified=False,
                       protected_sources=protected_sources)
        write_json(deployment_pending_path(context), pending)
    else:
        previous = pending.get("protected_sources")
        if (not isinstance(previous, dict) or previous.keys() != protected_sources.keys()
                or any(not persisted_identity_matches(previous[path], current)
                       for path, current in protected_sources.items())):
            raise ValueError("Pinned ClinVar source/index identity changed during recovery; preserved")
        record_retry_observation(pending, protected_sources=protected_sources)
        write_json(deployment_pending_path(context), pending)
    return retain_members(context, moves, reason=reason, phase="primary")


def retain_members(context: dict, moves: list, *, reason: str, phase: str) -> dict | None:
    """Journal the exact overwrite set before moving it; never delete bytes."""
    if not moves:
        return None
    root = context["root"].resolve()
    parent = root / ".installer-retained" / "clinvar-annotations"
    for directory in (parent.parent, parent):
        if directory.is_symlink() or (directory.exists() and not directory.is_dir()):
            raise ValueError(f"Unsafe ClinVar retention directory preserved: {directory}")
    parent.mkdir(parents=True, exist_ok=True)
    archive = Path(tempfile.mkdtemp(prefix=datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ-"), dir=parent))
    # Move the main annotation last; every overwrite member was preflighted.
    moves.sort(key=lambda item: item[0] == context["target"])
    journal = archive / "retirement.json"
    value = dict(kind="priva_clinvar_annotation_retirement_v1", status="retiring",
                 requested_assembly=context["assembly"], requested_vep_major=context["major"],
                 selected_vcf=str(context["target"]), reason=reason, created_utc=utc_now(),
                 scientific_qualified=False, phase=phase,
                 moves=[dict(source=str(path), retained=str(archive / path.name),
                             before=before, moved=False) for path, before in moves])
    write_json(journal, value)
    try:
        for move in value["moves"]:
            source, retained = Path(move["source"]), Path(move["retained"])
            if retirement_identity(source) != move["before"] or os.path.lexists(retained):
                raise ValueError("ClinVar retirement member changed or destination occupied; no rebuild started")
            source.rename(retained)
            after = retirement_identity(retained)
            move.update(moved=True, after=after)
            if any(after[key] != move["before"][key] for key in
                   ("dev", "ino", "mode", "nlink", "size", "mtime_ns")):
                raise ValueError("ClinVar member changed while being retained; no rebuild started")
            write_json(journal, value)
        value.update(status="retained_before_rebuild", completed_utc=utc_now())
        write_json(journal, value)
    except BaseException as error:
        value.update(status="retirement_incomplete", error=type(error).__name__ + ": " + str(error))
        write_json(journal, value)
        raise
    print("Retained previous managed ClinVar generation: " + str(journal), flush=True)
    return value


def retire_rebuilt_shards(config_path: Path, input_vcf: Path, chromosomes: list[str]) -> dict | None:
    """Internal hook after unchanged gnomAD query/selection, before shard reuse.

    Retrying after VEP or partial shard moves retains the newly built matching
    annotation and resumes only the remaining old shard paths.
    """
    context = deployment_context(config_path)
    pending = deployment_pending(context)
    if not context["owned"] or input_vcf.absolute() != context["target"]:
        raise ValueError("ClinVar recovery context does not match this selected input; no shards moved")
    if pending is None:
        return None
    before = retirement_identity(context["target"])
    deployment_header(context)
    if retirement_identity(context["target"]) != before:
        raise ValueError("New ClinVar VEP output changed before shard retirement; preserved")
    observed_sources = {}
    for path, recorded in pending["protected_sources"].items():
        current = retirement_identity(Path(path))
        if not persisted_identity_matches(recorded, current):
            raise ValueError("Pinned ClinVar source/index changed before shard retirement; preserved")
        observed_sources[path] = current
    if pending["phase"] in ("retiring_shards", "shards_retained"):
        if (not persisted_identity_matches(pending["rebuilt_target"], before)
                or pending["chromosomes"] != chromosomes):
            raise ValueError("New ClinVar VEP output/selected contigs changed during shard retirement; preserved")
        if pending["phase"] == "shards_retained":
            record_retry_observation(pending, protected_sources=observed_sources, rebuilt_target=before)
            write_json(deployment_pending_path(context), pending)
            return pending
    else:
        pending.update(phase="retiring_shards", rebuilt_target=before, chromosomes=chromosomes)
    members = shard_retirement_members(context, chromosomes)
    moves = retirement_inventory(context, members=members, protect_target=True)
    protected = [Path(path) for path in pending["protected_sources"]]
    if any(member.samefile(path) for member, _ in moves for path in protected):
        raise ValueError("A ClinVar shard aliases a supplied source/index; all shards preserved")
    record_retry_observation(pending, protected_sources=observed_sources, rebuilt_target=before)
    write_json(deployment_pending_path(context), pending)
    retain_members(context, moves, reason="cold VEP rebuild: retire mtime-only gnomAD reuse members", phase="shards")
    if retirement_identity(context["target"]) != before:
        raise ValueError("New ClinVar VEP output changed during shard retirement; preserved")
    pending.update(phase="shards_retained", shards_retained_utc=utc_now())
    write_json(deployment_pending_path(context), pending)
    return pending


def deploy_cache(config_path: Path, source_root: Path, *, source_vcf: Path | None = None,
                 runner=None) -> None:
    """Serialize one target independently of an enclosing gene-refresh lock."""
    from alphamissense_cache import AnnotationHeaderChangedError

    context = deployment_context(config_path)
    invoke = run_original if runner is None else runner
    arguments = [] if source_vcf is None else ["--source-vcf", source_vcf]
    if not context["owned"]:
        try:
            deployment_header(context)
        except (ValueError, EOFError, gzip.BadGzipFile, zlib.error) as error:
            raise ValueError(f"Explicit ClinVar annotation incompatible or unverifiable; preserved {context['target']}: {error}") from error
        invoke(source_root, "_ClinVar_VCF_deploy", config_path, *arguments)
        return
    # 2026-10-23 CT: never reuse the enclosing refresh workspace's .refresh.lock.
    lock_folder = context["folder"].resolve() / ".clinvar-deploy-locks" / context["target"].name
    lock_path = lock_folder / ".refresh.lock"
    if (lock_folder.parent.is_symlink() or lock_folder.is_symlink() or lock_path.is_symlink()
            or (lock_path.exists() and (not lock_path.is_file() or lock_path.stat().st_nlink != 1))):
        raise ValueError("Linked ClinVar deployment lock directory preserved")
    retirement_inventory(context, source_vcf)
    pinned_source_identity(context, source_vcf)
    with build_lock(lock_folder):
        # Waiting must not permit a different configuration/selection to inherit
        # the first target's lock or its managed ownership.
        if deployment_context(config_path) != context:
            raise ValueError("ClinVar deployment configuration changed while waiting for its lock")
        retirement_inventory(context, source_vcf)
        deployment_pending(context)
        if context["target"].exists():
            before = retirement_identity(context["target"])
            checked = {}
            try:
                deployment_header(context)
                # The shared validity helper can repair a damaged gzip from a
                # newer plain sibling. Check that possible source beforehand,
                # retaining both generations if its provenance is incompatible.
                plain = Path(str(context["target"])[:-3])
                if plain.exists() and plain.stat().st_mtime_ns > context["target"].stat().st_mtime_ns:
                    checked[plain] = retirement_identity(plain)
                    deployment_header(context, plain)
            except AnnotationHeaderChangedError:
                raise
            except (ValueError, EOFError, gzip.BadGzipFile, zlib.error) as error:
                retire_owned_cache(config_path, reason=str(error), source_vcf=source_vcf,
                                   expected_target=before, expected_members=checked)
        invoke(source_root, "_ClinVar_VCF_deploy", config_path, *arguments)
        deployment_header(deployment_context(config_path))
        pending = deployment_pending_path(context)
        if pending.exists():
            pending.unlink()


def require_workspace_owner(plan: dict) -> None:
    """An interrupted annotation is reusable only within its recorded profile."""
    folder = Path(plan["build_dir"])
    previous_path = folder / "refresh.plan.json"
    if previous_path.exists():
        previous = json.loads(previous_path.read_text())
        if (previous.get("kind") != RECEIPT_KIND
                or profile_signature(previous.get("scientific_profile")) != profile_signature(plan["scientific_profile"])
                or previous.get("explicit_vcf_identity") != plan["explicit_vcf_identity"]
                or previous.get("latest", {}).get("file_date") != plan["latest"]["file_date"]
                or previous.get("latest", {}).get("publisher_md5") != plan["latest"]["publisher_md5"]):
            raise ValueError("Occupied incomplete ClinVar workspace belongs to a different or unknown scientific profile/input; "
                             "its annotation and previous plan were retained without reuse: " + str(folder))
    else:
        occupied = ((folder / "install.config.yaml").exists()
                    or Path(plan["output_path"]).exists()
                    or any(folder.glob("clinvar.*.vep*")))
        if occupied:
            raise ValueError("Occupied ClinVar annotation workspace has no owned scientific-profile plan; "
                             "it was retained without reuse: " + str(folder))


def complete_refresh(plan: dict, working_config: Path, inputs: dict) -> None:
    original = Path(plan["config_path"])
    if digest(original) != plan["config_sha256"]:
        raise ValueError("The original configuration changed during ClinVar regeneration; no new selection was written")
    baseline = Path(plan["baseline"]["path"])
    before = {key: value for key, value in plan["baseline"].items() if key != "refresh_boundary_utc_date"}
    if baseline_metadata(baseline) != before:
        raise ValueError("Bundled ClinVar metadata changed during regeneration; no new selection was written")
    current = load_config(working_config)
    current_policy = resource_policy(Path(plan["source_root"]))
    if profile_signature(scientific_profile(load_config(original), Path(plan["source_root"]), current_policy)) != profile_signature(plan["scientific_profile"]):
        raise ValueError("ClinVar scientific inputs/source changed during regeneration; no new selection was written")
    vcf = selected_vcf(current, Path(current["clinvar_vcf_dir"]), plan["ucsc_assembly"])
    if local_file_date(vcf).isoformat() != plan["latest"]["file_date"]:
        raise ValueError("Completed annotation fileDate differs from the detected latest ClinVar release")
    updates = {key: str(current[key]) for key in ("clinvar_vcf", "clinvar_aa_stat", "clinvar_splice_stat")}
    updates["clinvar_gene_stat"] = plan["output_path"]
    output = Path(plan["output_path"])
    value = {"kind": RECEIPT_KIND, "status": "complete", "data_date": plan["latest"]["file_date"],
             "scientific_qualified": False,
             "assembly": plan["assembly"], "vep_major": plan["vep_major"], "baseline": plan["baseline"],
             "managed_root": plan["managed_root"], "latest": plan["latest"], "inputs": inputs,
             "resource_policy": plan["resource_policy"],
             "scientific_profile": plan["scientific_profile"],
             "scientific_profile_sha256": plan["scientific_profile_sha256"],
             "completed_artifacts": {"vcf": identity(vcf), "vcf_index": identity(Path(str(vcf) + ".tbi")),
                                     "aa": identity(Path(current["clinvar_aa_stat"])),
                                     "splice": identity(Path(current["clinvar_splice_stat"]))},
             "output": identity(output), "annotation": identity(vcf), "completed_utc": utc_now(),
             "config_updates": updates, "source_root": plan["source_root"],
             "scientific_sources": {name: digest(Path(plan["source_root"]) / "scripts" / name)
                                    for name in ("clinvar_stat_variants.py", "combine_annotations.py", "annotation_vcf.sh")},
             "qualification": "software refresh; public/family acceptance remains separate"}
    write_json(receipt_path(output), value)
    atomic_text(original, render_config(original.read_text(), updates))


def install_refresh(config_path: Path, source_root: Path) -> dict:
    # Preflight selected inputs before metadata access or any acquisition.
    validate_main(load_config(config_path), source_root, allow_missing_bundled_clinvar=True)
    plan = plan_refresh(config_path, source_root)
    print("ClinVar gene-statistics refresh: " + plan["action"] + "; " + plan["reason"], flush=True)
    if plan["action"] == "preserve":
        for function in ("ClinVar_VCF_deploy", "ClinVar_AA_stat", "ClinVar_Gene_stat"):
            run_original(source_root, function, config_path)
        return plan
    folder = Path(plan["build_dir"])
    with build_lock(folder):
        prior = owned_receipt(Path(plan["output_path"]), assembly=plan["assembly"], vep_major=plan["vep_major"])
        if prior is not None and prior["data_date"] == plan["latest"]["file_date"]:
            validate_reuse(prior, plan)
            if digest(config_path) != plan["config_sha256"]:
                raise ValueError("Configuration changed while waiting for the ClinVar refresh lock")
            atomic_text(config_path, render_config(config_path.read_text(), prior["config_updates"]))
            return dict(plan, action="reuse")
        require_workspace_owner(plan)
        write_json(folder / "refresh.plan.json", plan)
        inputs = acquire_inputs(plan)
        working = folder / "install.config.yaml"
        updates = {}
        if not plan["explicit_vcf"]:
            updates = {"clinvar_vcf_dir": str(folder),
                       "clinvar_vcf": str(folder / f"clinvar.{plan['ucsc_assembly']}.vep.vcf.gz")}
            for key, suffix in (("clinvar_aa_stat", "aa_change.stats.pkl"),
                                ("clinvar_splice_stat", "splice_change.stats.pkl")):
                updates[key] = str(folder / f"clinvar.{plan['ucsc_assembly']}.vep.{suffix}")
        atomic_text(working, render_config(config_path.read_text(), updates))
        arguments = [] if plan["explicit_vcf"] else ["--source-vcf", inputs["sources"]["vcf"]["path"]]
        run_original(source_root, "ClinVar_VCF_deploy", working, *arguments)
        annotated = selected_vcf(load_config(working), Path(load_config(working)["clinvar_vcf_dir"]), plan["ucsc_assembly"])
        if local_file_date(annotated).isoformat() != plan["latest"]["file_date"]:
            raise ValueError("Managed ClinVar annotation has an older/different fileDate; refusing reuse")
        run_original(source_root, "ClinVar_AA_stat", working)
        run_original(source_root, "ClinVar_Gene_stat_regenerate", working,
                     inputs["sources"]["summary"]["path"], inputs["sources"]["mane"]["path"], plan["output_path"])
        complete_refresh(plan, working, inputs)
        return plan


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    commands = parser.add_subparsers(dest="command", required=True)
    for name in ("plan", "install", "deploy", "retire-owned", "retire-shards", "check-header"):
        command = commands.add_parser(name)
        command.add_argument("--config", type=Path, required=True)
        command.add_argument("--source-root", type=Path, required=True)
        if name == "plan":
            command.add_argument("--output", type=Path)
        if name in {"deploy", "retire-owned"}:
            command.add_argument("--source-vcf", type=Path)
        if name == "retire-shards":
            command.add_argument("--input-vcf", type=Path, required=True)
            command.add_argument("--chromosome", action="append", required=True)
    command = commands.add_parser("acquire")
    command.add_argument("--plan", type=Path, required=True)
    command.add_argument("--output", type=Path, required=True)
    command = commands.add_parser("file-date")
    command.add_argument("--input", type=Path, required=True)
    command = commands.add_parser("candidate-route", help="Check preparation source wiring without reading cache inputs")
    command.add_argument("--source-root", type=Path, required=True)
    for name in ("candidate-begin", "candidate-complete"):
        command = commands.add_parser(name, help="Record an explicitly unqualified managed ClinVar VEP preparation")
        command.add_argument("--config", type=Path, required=True)
        command.add_argument("--source-root", type=Path, required=True)
        command.add_argument("--input-vcf", type=Path, required=True)
    args = parser.parse_args()
    try:
        if args.command == "deploy":
            deploy_cache(args.config, args.source_root, source_vcf=args.source_vcf)
        elif args.command == "retire-owned":
            retire_owned_cache(args.config, reason="existing deployment requires a cold annotation rebuild",
                               source_vcf=args.source_vcf)
        elif args.command == "retire-shards":
            retire_rebuilt_shards(args.config, args.input_vcf, args.chromosome)
        elif args.command == "check-header":
            print(json.dumps(deployment_header(deployment_context(args.config)), sort_keys=True))
        elif args.command == "candidate-route":
            print(json.dumps(candidate_route(args.source_root), sort_keys=True))
        elif args.command in ("candidate-begin", "candidate-complete"):
            function = candidate_begin if args.command == "candidate-begin" else candidate_complete
            print(json.dumps(function(args.config, args.source_root, args.input_vcf), sort_keys=True))
        elif args.command == "file-date":
            print(local_file_date(args.input).isoformat())
        elif args.command == "acquire":
            write_json(args.output, acquire_inputs(json.loads(args.plan.read_text())))
        elif args.command == "plan":
            value = plan_refresh(args.config, args.source_root)
            if args.output:
                write_json(args.output, value)
            else:
                print(json.dumps(value, indent=2, sort_keys=True))
        else:
            install_refresh(args.config, args.source_root)
    except (OSError, EOFError, ValueError, KeyError, zlib.error, subprocess.CalledProcessError) as error:
        parser.exit(1, f"ERROR: ClinVar gene-statistics refresh failed: {error}\n")


if __name__ == "__main__":
    main()
