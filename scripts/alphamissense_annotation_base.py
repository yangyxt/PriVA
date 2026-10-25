#!/usr/bin/env python3
"""Keep the complete AlphaMissense VEP input separate from its gnomAD derivative.

The scientific converter and basic_vep_annotation procedure remain in the shell
installer. This helper records their real inputs and publishes completed files;
it never changes their flags, row selection, scoring or mutation order.

2026-10-22 CT: the original installer computes motifs before gnomAD annotation.
The latter intentionally retains only its supported chromosomes, including no
chrM. A timestamp on that derivative cannot identify the complete motif input.
"""
from __future__ import annotations

import argparse
from contextlib import contextmanager
import fcntl
import hashlib
import importlib.util
import importlib
import importlib.metadata
import json
import os
from pathlib import Path
import re
import shutil
import stat
import subprocess
import sys
import tempfile

import yaml

SPEC = importlib.util.spec_from_file_location("alpha_base_header", Path(__file__).with_name("alphamissense_cache.py"))
HEADER = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(HEADER)
SCHEMA = "priva-alphamissense-annotation-base-v1"
MOTIF_SCHEMA = "priva-alphamissense-motif-input-v1"


def require(condition, message):
    if not condition:
        raise ValueError(message)


def digest(path):
    result = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for block in iter(lambda: stream.read(8 * 1024 * 1024), b""):
            result.update(block)
    return result.hexdigest()


def json_digest(value):
    return hashlib.sha256(json.dumps(value, sort_keys=True, separators=(",", ":")).encode()).hexdigest()


def identity(path, missing=False):
    path = Path(path).absolute()
    if missing and not os.path.lexists(path):
        return None
    info = path.stat()
    require(stat.S_ISREG(info.st_mode) and (missing or info.st_size > 0), "Missing/empty/nonregular dependency: " + str(path))
    # st_dev varies between cluster mounts of the same native filesystem.
    return dict(path=str(path), target=str(path.resolve()), inode=info.st_ino,
                size=info.st_size, mtime_ns=info.st_mtime_ns, ctime_ns=info.st_ctime_ns)


def tree_identity(directory, code=False):
    """Metadata for the selected cache; content hashes for the actual VEP API.

    No unrelated VEP release, plugin, or whole Conda environment is scanned.
    Cache payload bytes are not reread on a warm check.
    """
    directory = Path(directory).resolve()
    require(directory.is_dir(), "Missing selected VEP dependency tree: " + str(directory))
    result = []
    total_bytes = latest_mtime = latest_ctime = 0
    for root, names, files in os.walk(directory):
        names[:] = sorted(name for name in names if name not in (".git", "__pycache__"))
        for name in sorted(files):
            path = Path(root) / name
            info = path.stat()
            require(stat.S_ISREG(info.st_mode), "Nonregular selected VEP dependency: " + str(path))
            item = dict(name=str(path.relative_to(directory)),
                        size=info.st_size, inode=info.st_ino,
                        mtime_ns=info.st_mtime_ns, ctime_ns=info.st_ctime_ns)
            if code:
                # API source is small; its content, rather than installation
                # timestamps, permits equivalent original/Nextflow checkouts.
                item = dict(name=item["name"], sha256=digest(path))
            result.append(item)
            total_bytes += info.st_size
            latest_mtime = max(latest_mtime, info.st_mtime_ns)
            latest_ctime = max(latest_ctime, info.st_ctime_ns)
    require(result, "Empty selected VEP dependency tree: " + str(directory))
    return dict(path=str(directory), files=len(result), bytes=total_bytes,
                inventory_sha256=json_digest(result), latest_mtime_ns=latest_mtime, latest_ctime_ns=latest_ctime)


def procedure_text(path, name):
    text = Path(path).read_text()
    match = re.search(r"^function " + re.escape(name) + r"\(\) \{\n.*?^\}", text, re.M | re.S)
    require(match is not None, "Cannot identify scientific shell procedure: " + name)
    return match.group()


def shell_procedure(path, name):
    return hashlib.sha256(procedure_text(path, name).encode()).hexdigest()


def annotation_recipe(path):
    body = procedure_text(path, "AlphaMissense_anno")
    converter = re.search(r'^\s*python "\$convert_py"[^\n]+', body, re.M)
    call = re.search(r'^[ \t]*basic_vep_annotation \\\n', body, re.M)
    require(converter is not None and call is not None, "Cannot identify actual AlphaMissense converter/VEP invocation")
    lines = []
    for line in body[call.start():].splitlines():
        lines.append(line)
        if re.search(r'\s&&\s*\\?$', line) or not line.rstrip().endswith("\\"):
            break
    return hashlib.sha256((converter.group().strip() + "\n" + "\n".join(lines).strip()).encode()).hexdigest()


def vep_execution(path):
    body = procedure_text(path, "basic_vep_annotation")
    command = re.search(r'^\s*vep\s+[^\n]+', body, re.M)
    require(command is not None, "Cannot identify the actual VEP execution block")
    plugin = body.find('    if [[ -f ${spliceai_snv_prescore}')
    start = plugin if plugin >= 0 else command.start()
    # Config-file quoting and log text are provenance; actual optional plugin
    # construction, every VEP flag and sort/index command remain scientific.
    execution = "\n".join(line for line in body[start:].splitlines() if not line.lstrip().startswith("log "))
    return hashlib.sha256(execution.encode()).hexdigest()


def settings(config_file, source_root, base=None):
    config = yaml.safe_load(Path(config_file).read_text())
    require(isinstance(config, dict), "Installer configuration must be a mapping")
    root = Path(config.get("cache_root") or "")
    require(root.is_absolute(), "Configure the managed cache_root before AlphaMissense annotation")
    canonical = Path(config["alphamissense_vep_vcf"])
    base = Path(base or config.get("alphamissense_vep_base_vcf") or
                canonical.with_name(canonical.name.removesuffix(".vcf.gz") + ".base.vcf.gz"))
    paths = dict(converted=Path(config["alphamissense_vcf"]), base=base, canonical=canonical)
    for path in paths.values():
        require(path.is_absolute() and path.name.endswith(".vcf.gz") and ".." not in path.parts,
                "Use absolute managed AlphaMissense .vcf.gz destinations")
        require(path.resolve().is_relative_to(root.resolve()) and not path.is_symlink(),
                "AlphaMissense destinations must be independent files inside cache_root: " + str(path))
    selected = [path.resolve() for path in paths.values()]
    require(len(set(selected)) == len(selected), "Converted, complete base and runtime VCF paths must differ")
    present = [path for path in paths.values() if path.exists()]
    require(all(not left.samefile(right) for i, left in enumerate(present) for right in present[i + 1:]),
            "AlphaMissense base and runtime files must not be hardlinks")
    return config, Path(source_root), paths


def science_inputs(config, source_root):
    assembly = HEADER.assembly_name(config["assembly"])
    genome = "GRCh37" if assembly == "hg19" else "GRCh38"
    major = str(config["vep_installed_version"])
    require(re.fullmatch(r"[1-9][0-9]*", major), "VEP major must be numeric")
    executable = shutil.which("vep")
    require(executable, "Activate the configured PriVA environment before recording annotation inputs")
    executable = Path(executable).resolve()
    help_text = subprocess.run([str(executable), "--help"], check=True, text=True,
                               stdout=subprocess.PIPE, stderr=subprocess.STDOUT).stdout
    version = re.search(r"ensembl-vep\s*:?\s+([0-9]+(?:\.[0-9]+)*)", help_text)
    require(version is not None and version.group(1).split(".")[0] == major,
            "Runnable VEP does not match the configured major")
    api = executable.parent
    # Match the executable's actual FindBin/lib search order, then ask Perl
    # where it loaded the modules. PERL5LIB can otherwise select a different API
    # even when the executable and its apparent sibling directories look right.
    perl = shutil.which("perl")
    require(perl, "The active VEP Perl runtime is missing")
    conservation = config.get("conservation_file")
    plugin = Path(config["vep_plugins_dir"]) / "Conservation.pm"
    query = ('if (@ARGV) { require $ARGV[0]; } my %loaded = %INC; '
        'my @shared = @DynaLoader::dl_shared_objects; require JSON::PP; '
        'print JSON::PP::encode_json({version => $], shared => [@shared], '
        'loaded => {map {$_ => $loaded{$_}} grep {defined $loaded{$_} && !ref $loaded{$_}} keys %loaded}})')
    loaded = json.loads(subprocess.run([perl, "-I" + str(api), "-I" + str(api / "modules"),
        "-MGetopt::Long", "-MFindBin", "-MBio::EnsEMBL::VEP::Runner", "-MBio::EnsEMBL::VEP::Utils", "-e", query] +
        ([] if conservation in (None, "", "null") else [str(plugin)]),
        check=True, text=True, stdout=subprocess.PIPE).stdout)
    require("Bio/EnsEMBL/VEP/Runner.pm" in loaded["loaded"], "Cannot observe the API actually loaded by VEP")
    api_roots = set()
    for name, location in loaded["loaded"].items():
        if name.startswith("Bio/EnsEMBL/"):
            path = Path(location).resolve()
            for _ in Path(name).parts:
                path = path.parent
            api_roots.add(path / "Bio")
    api_trees = [tree_identity(path, code=True) for path in sorted(api_roots)]
    loaded_files = {name: dict(path=str(Path(path).resolve()), sha256=digest(path))
                    for name, path in sorted(loaded["loaded"].items()) if Path(path).is_file()}
    for path in loaded.get("shared", []):
        loaded_files["XS:" + str(Path(path).resolve())] = dict(path=str(Path(path).resolve()), sha256=digest(path))
    cache = Path(config["vep_cache_dir"]) / "homo_sapiens_merged" / f"{major}_{genome}"
    reference = Path(config["ref_genome"])
    return dict(assembly=assembly, vep_version=version.group(1), threads=int(config["threads"]),
                converter_sha256=digest(source_root / "scripts/alphmis_tsv2vcf.py"),
                basic_vep_execution_sha256=vep_execution(source_root / "scripts/install_utils.sh"),
                annotation_recipe_sha256=annotation_recipe(source_root / "scripts/install_utils.sh"),
                prescore=[identity(config["alphamissense_prescore"]), identity(str(config["alphamissense_prescore"]) + ".tbi")],
                reference=[identity(reference), identity(str(reference) + ".fai")] +
                    ([identity(str(reference) + ".gzi")] if Path(str(reference) + ".gzi").exists() else []),
                conservation=None if conservation in (None, "", "null") else
                    dict(track=identity(conservation), plugin_sha256=digest(plugin)),
                vep=dict(executable=str(executable), sha256=digest(executable), api=api_trees,
                         perl=identity(perl), perl_version=loaded["version"],
                         loaded=dict(files=len(loaded_files), inventory_sha256=json_digest(loaded_files),
                                     runner_path=str(Path(loaded["loaded"]["Bio/EnsEMBL/VEP/Runner.pm"]).resolve()))),
                cache=tree_identity(cache))


def receipt_path(base):
    return Path(str(base) + ".annotation-base.json")


def pending_path(base):
    return Path(str(base) + ".annotation-base.pending.json")


def output_state(paths):
    return {key: identity(path, missing=True) for key, path in
            (("converted", paths["converted"]), ("converted_index", Path(str(paths["converted"]) + ".tbi")),
             ("base", paths["base"]), ("base_index", Path(str(paths["base"]) + ".tbi")))}


def atomic_json(path, value):
    fd, temporary = tempfile.mkstemp(prefix="." + path.name + ".", dir=path.parent)
    try:
        with os.fdopen(fd, "w") as stream:
            json.dump(value, stream, sort_keys=True, indent=2)
            stream.write("\n")
            stream.flush()
            os.fsync(stream.fileno())
        os.replace(temporary, path)
    finally:
        if os.path.exists(temporary):
            os.unlink(temporary)


def index_counts(path):
    result = subprocess.run(["bcftools", "index", "--stats", str(path)], check=True,
                            text=True, stdout=subprocess.PIPE).stdout
    rows = [line.split("\t") for line in result.splitlines()]
    require(rows and all(len(row) == 3 and row[2].isdigit() for row in rows),
            "Missing complete per-contig VCF index counts: " + str(path))
    return {row[0]: int(row[2]) for row in rows if int(row[2])}


@contextmanager
def directory_lock(directory):
    descriptor = os.open(directory, os.O_RDONLY | os.O_DIRECTORY)
    try:
        fcntl.flock(descriptor, fcntl.LOCK_EX | fcntl.LOCK_NB)
        yield
    finally:
        os.close(descriptor)


def capture(config_file, source_root, base=None):
    config, source_root, paths = settings(config_file, source_root, base)
    return dict(science=science_inputs(config, source_root), prior=output_state(paths),
                runtime_prior=[identity(paths["canonical"], missing=True), identity(str(paths["canonical"]) + ".tbi", missing=True)])


def check(config_file, source_root, base=None, pending=False):
    config, source_root, paths = settings(config_file, source_root, base)
    require(pending or not pending_path(paths["base"]).exists(), "Complete the pending base publication first")
    record = json.loads(receipt_path(paths["base"]).read_text())
    require(record.get("schema") == SCHEMA and record.get("science") == science_inputs(config, source_root),
            "Complete annotation producer or scientific dependency changed")
    require(record.get("outputs") == output_state(paths), "Complete annotation files changed or publication is incomplete")
    header = HEADER.check_annotation(input_file=paths["base"], assembly=config["assembly"], vep_major=config["vep_installed_version"])
    require(header["vep_version"] == record["science"]["vep_version"], "Complete annotation VEP patch differs")
    require(isinstance(record.get("base_sha256"), str) and re.fullmatch(r"[0-9a-f]{64}", record["base_sha256"]),
            "Complete annotation has no published content identity")
    return record


def moved_identity(actual, expected):
    return actual is not None and expected is not None and all(
        actual[key] == expected[key] for key in ("inode", "size", "mtime_ns"))


def commit(config_file, source_root, base=None):
    config, source_root, paths = settings(config_file, source_root, base)
    marker = pending_path(paths["base"])
    require(marker.is_file() and not marker.is_symlink() and marker.stat().st_size < 4 * 1024 * 1024,
            "Missing or invalid base publication marker")
    value = json.loads(marker.read_text())
    require(value.get("schema") == SCHEMA and value.get("science") == science_inputs(config, source_root),
            "Pending annotation scientific inputs changed")
    require(value.get("targets") == {key: str(path) for key, path in paths.items()}, "Pending publication paths changed")
    if value.get("native_runtime") is not None:
        validate_runtime_lineage(value["native_runtime"], paths, value["science"], value["base_sha256"])
        require(value["preserved_converted"] == [identity(paths["converted"]), identity(str(paths["converted"]) + ".tbi")],
                "Preserved native converted files changed during import")
    locations = []
    root = Path(config["cache_root"]).resolve()
    for operation in value["operations"]:
        source, target = Path(operation["source"]), Path(operation["target"])
        require(source.resolve().is_relative_to(root) and target.resolve().is_relative_to(root)
                and not source.is_symlink() and not target.is_symlink(), "Pending annotation path escapes managed storage")
        if source.exists():
            require(identity(source) == operation["new"] and identity(target, missing=True) == operation["prior"],
                    "Staged annotation or previous target changed")
            locations.append("stage")
        else:
            require(moved_identity(identity(target, missing=True), operation["new"]), "Published annotation changed")
            locations.append("final")
    require(locations == sorted(locations), "Pending publication is not an ordered prefix")
    for location, operation in zip(locations, value["operations"]):
        if location == "stage":
            os.replace(operation["source"], operation["target"])
    if value["kind"] == "base":
        atomic_json(receipt_path(paths["base"]), dict(schema=SCHEMA, science=value["science"],
                    outputs=output_state(paths), base_sha256=value["base_sha256"], counts=value["counts"],
                    lineage=value.get("lineage"), native_runtime=value.get("native_runtime"),
                    producer_paths=value.get("producer_paths")))
    marker.unlink()


def recover(config_file, source_root, base=None):
    _, _, paths = settings(config_file, source_root, base)
    if pending_path(paths["base"]).exists():
        with directory_lock(paths["base"].parent):
            commit(config_file, source_root, base)


def copy_verified(source, target, expected=None):
    before = identity(source)
    result = hashlib.sha256()
    with Path(source).open("rb") as original, Path(target).open("xb") as output:
        for block in iter(lambda: original.read(8 * 1024 * 1024), b""):
            output.write(block)
            result.update(block)
        output.flush()
        os.fsync(output.fileno())
    require(identity(source) == before, "Annotation source changed while copying")
    value = result.hexdigest()
    require(expected is None or value == expected, "Verified native source SHA256 differs")
    return value


def publish(config_file, source_root, staged_converted, staged_annotation, expected_inputs, base=None, lineage=None):
    config, source_root, paths = settings(config_file, source_root, base)
    expected = json.loads(expected_inputs.read_text()) if isinstance(expected_inputs, Path) else (
        json.loads(expected_inputs) if isinstance(expected_inputs, str) else expected_inputs)
    require(capture(config_file, source_root, base) == expected, "Annotation inputs or previous outputs changed during generation")
    converted, annotation = Path(staged_converted), Path(staged_annotation)
    require(converted != annotation and all(path.resolve().is_relative_to(Path(config["cache_root"]).resolve())
            for path in (converted, annotation)), "Scientific stages must remain inside cache_root")
    counts = index_counts(converted)
    require(counts == index_counts(annotation), "VEP output lost or added input records, including unsupported AF chromosomes")
    header = HEADER.check_annotation(input_file=annotation, assembly=config["assembly"], vep_major=config["vep_installed_version"])
    require(header["vep_version"] == expected["science"]["vep_version"], "Generated VEP patch differs from actual producer")
    for path in paths.values():
        path.parent.mkdir(parents=True, exist_ok=True)
    with directory_lock(paths["base"].parent):
        require(not pending_path(paths["base"]).exists(), "Recover the pending base publication before rebuilding")
        require(capture(config_file, source_root, base) == expected, "Annotation inputs changed before publication")
        stage = Path(tempfile.mkdtemp(prefix=".alpha-base-runtime-", dir=paths["canonical"].parent))
        runtime = stage / paths["canonical"].name
        base_sha = copy_verified(annotation, runtime)
        copy_verified(str(annotation) + ".tbi", str(runtime) + ".tbi")
        operations = []
        for source, target in ((converted, paths["converted"]), (annotation, paths["base"]), (runtime, paths["canonical"])):
            for suffix in ("", ".tbi"):
                source_file, target_file = Path(str(source) + suffix), Path(str(target) + suffix)
                require(source_file != target_file, "Publication needs an independent prepared source")
                operations.append(dict(source=str(source_file), target=str(target_file), new=identity(source_file),
                                       prior=identity(target_file, missing=True)))
        atomic_json(pending_path(paths["base"]), dict(schema=SCHEMA, kind="base", science=expected["science"],
                    targets={key: str(path) for key, path in paths.items()}, operations=operations,
                    base_sha256=base_sha, counts=counts, lineage=lineage,
                    producer_paths=dict(source_root=str(source_root.absolute()), converter=str(source_root / "scripts/alphmis_tsv2vcf.py"),
                                        basic_vep_annotation=str(source_root / "scripts/install_utils.sh"),
                                        full_basic_procedure_sha256=shell_procedure(source_root / "scripts/install_utils.sh", "basic_vep_annotation"))))
        commit(config_file, source_root, base)


def restore(config_file, source_root, base=None):
    record = check(config_file, source_root, base)
    _, _, paths = settings(config_file, source_root, base)
    paths["canonical"].parent.mkdir(parents=True, exist_ok=True)
    with directory_lock(paths["base"].parent):
        stage = Path(tempfile.mkdtemp(prefix=".alpha-base-runtime-", dir=paths["canonical"].parent))
        runtime = stage / paths["canonical"].name
        copy_verified(paths["base"], runtime, record["base_sha256"])
        copy_verified(str(paths["base"]) + ".tbi", str(runtime) + ".tbi")
        operations = [dict(source=str(source), target=str(target), new=identity(source), prior=identity(target, missing=True))
                      for source, target in ((runtime, paths["canonical"]),
                          (Path(str(runtime) + ".tbi"), Path(str(paths["canonical"]) + ".tbi")))]
        atomic_json(pending_path(paths["base"]), dict(schema=SCHEMA, kind="runtime", science=record["science"],
                    targets={key: str(path) for key, path in paths.items()}, operations=operations))
        commit(config_file, source_root, base)


def adopt(config_file, source_root, lineage_file, base=None):
    """Explicit migration of externally verified *complete* native artifacts.

    A lineage document supplies schema, current captured science, per-contig
    counts, provenance, and converted/annotation {path,sha256,index_sha256}.
    The caller must obtain that lineage from the actual completed native run;
    a recent/enriched VCF is never automatically promoted into a base.
    """
    lineage = json.loads(Path(lineage_file).read_text())
    expected = capture(config_file, source_root, base)
    require(lineage.get("schema") == SCHEMA and lineage.get("science") == expected["science"]
            and lineage.get("provenance"), "Native complete-base lineage does not match the actual scientific producer")
    _, _, paths = settings(config_file, source_root, base)
    paths["base"].parent.mkdir(parents=True, exist_ok=True)
    stage = Path(tempfile.mkdtemp(prefix=".alpha-base-adopt-", dir=paths["base"].parent))
    native_runtime = lineage.get("native_runtime")
    if native_runtime is not None:
        validate_runtime_lineage(native_runtime, paths, expected["science"], lineage["annotation"]["sha256"])
        # Closed native import adds only the new full base, preserving the
        # observed converted/runtime bytes and timestamps. The runtime may be
        # pre-AF or AF-enriched; normal installation still completes unfinished
        # AF work. This import does not imply whole-installer success.
        converted = Path(lineage["converted"]["path"])
        require(converted.resolve() == paths["converted"].resolve(), "Preserved native converted target differs")
        before = identity(converted), identity(str(converted) + ".tbi")
        require(digest(converted) == lineage["converted"]["sha256"] and
                digest(str(converted) + ".tbi") == lineage["converted"]["index_sha256"] and
                index_counts(converted) == lineage["counts"], "Preserved native converted identity/counts differ")
        annotation = stage / "annotation.vcf.gz"
        source = lineage["annotation"]
        base_sha = copy_verified(source["path"], annotation, source["sha256"])
        copy_verified(source["path"] + ".tbi", str(annotation) + ".tbi", source["index_sha256"])
        require(index_counts(annotation) == lineage["counts"], "Native complete-base chromosome counts differ")
        header = HEADER.check_annotation(input_file=annotation, assembly=expected["science"]["assembly"],
                                         vep_major=expected["science"]["vep_version"].split(".")[0])
        require(header["vep_version"] == expected["science"]["vep_version"], "Native complete-base VEP patch differs")
        with directory_lock(paths["base"].parent):
            require(capture(config_file, source_root, base) == expected and
                    before == (identity(converted), identity(str(converted) + ".tbi")), "Native import inputs changed")
            validate_runtime_lineage(native_runtime, paths, expected["science"], base_sha)
            require(not pending_path(paths["base"]).exists(), "Recover the pending base import first")
            operations = [dict(source=str(src), target=str(dst), new=identity(src), prior=identity(dst, missing=True))
                for src, dst in ((annotation, paths["base"]), (Path(str(annotation) + ".tbi"), Path(str(paths["base"]) + ".tbi")))]
            atomic_json(pending_path(paths["base"]), dict(schema=SCHEMA, kind="base", science=expected["science"],
                targets={key: str(path) for key, path in paths.items()}, operations=operations, base_sha256=base_sha,
                counts=lineage["counts"], native_runtime=native_runtime,
                preserved_converted=list(before),
                producer_paths=dict(source_root=str(source_root.absolute())),
                lineage=dict(path=str(Path(lineage_file).absolute()), sha256=digest(lineage_file))))
            commit(config_file, source_root, base)
        return
    selected = []
    for key in ("converted", "annotation"):
        source = lineage[key]
        target = stage / (key + ".vcf.gz")
        copy_verified(source["path"], target, source["sha256"])
        copy_verified(source["path"] + ".tbi", str(target) + ".tbi", source["index_sha256"])
        require(index_counts(target) == lineage["counts"], "Native complete-base chromosome counts differ")
        selected.append(target)
    publish(config_file, source_root, *selected, expected, base, dict(path=str(Path(lineage_file).absolute()), sha256=digest(lineage_file)))


def validate_runtime_lineage(value, paths, science, base_sha):
    require(isinstance(value, dict) and value.get("vcf") == identity(paths["canonical"])
            and value.get("index") == identity(str(paths["canonical"]) + ".tbi")
            and value.get("annotation_base_sha256") == base_sha
            and value.get("science_sha256") == json_digest(science),
            "Closed native runtime is not bound to this exact complete annotation and scientific profile")
    proof = value.get("closed_native_proof", {})
    require(proof.get("path") and proof.get("sha256") == digest(proof["path"]),
            "Closed native runtime proof is missing or changed")


def runtime_associated(config_file, source_root, base=None):
    record = check(config_file, source_root, base)
    _, _, paths = settings(config_file, source_root, base)
    validate_runtime_lineage(record.get("native_runtime"), paths, record["science"], record["base_sha256"])


def motif_binding(config_file, source_root, base=None):
    record = check(config_file, source_root, base)
    return dict(schema=MOTIF_SCHEMA, annotation_base_sha256=record["base_sha256"],
                annotation_science_sha256=json_digest(record["science"]), counts=record["counts"],
                producer_sha256=digest(Path(source_root) / "scripts/am_pick_missense_constrained_segments.py"),
                runtime=motif_runtime())


def motif_runtime():
    """Only the unchanged scientific producer's Python and imported packages.

    File metadata catches package replacement without rehashing unrelated
    environments. Version/path inventories and latest ctime are also retained
    so an explicit historical import can state what is actually evidenced;
    this function cannot invent a package snapshot for an older native run.
    """
    packages = {}
    for module_name, distribution_name in (("numpy", "numpy"), ("pandas", "pandas"), ("scipy", "scipy"),
            ("sklearn", "scikit-learn"), ("pysam", "pysam"), ("joblib", "joblib")):
        module = importlib.import_module(module_name)
        distribution = importlib.metadata.distribution(distribution_name)
        owned = []
        for item in distribution.files or ():
            if "__pycache__" in item.parts or str(item).endswith((".pyc", ".pyo")):
                continue
            path = Path(distribution.locate_file(item)).absolute()
            try:
                value = path.stat()
            except FileNotFoundError:
                continue
            if stat.S_ISREG(value.st_mode):
                owned.append(dict(name=str(item), path=str(path), inode=value.st_ino, size=value.st_size,
                                  mtime_ns=value.st_mtime_ns, ctime_ns=value.st_ctime_ns))
        require(owned, "Cannot identify scientific Python package files: " + distribution_name)
        packages[module_name] = dict(version=str(module.__version__), distribution_version=distribution.version,
            module_path=str(Path(module.__file__).resolve()), files=len(owned),
            inventory_sha256=json_digest(sorted(owned, key=lambda value: value["name"])),
            latest_ctime_ns=max(value["ctime_ns"] for value in owned))
    return dict(python=identity(sys.executable), version=sys.version, packages=packages)


def motif_check(config_file, source_root, artifact, base=None):
    record = json.loads(Path(str(artifact) + ".nm_to_enst.json").read_text())
    require(record.get("scientific_input") == motif_binding(config_file, source_root, base),
            "Motif receipt is not bound to this complete annotation and scientific producer")


def adopt_motif(config_file, source_root, artifact, lineage_file, mapping, inputs, base=None):
    """Bind a preserved, natively traced first motif without running CV again.

    The explicit closed-run lineage must pin its full base SHA, scientific
    input/producer binding and exact artifact/mapping identity. Mtime or an
    identical-looking current output cannot supply that evidence.
    """
    from types import SimpleNamespace
    spec = importlib.util.spec_from_file_location("alpha_motif_publisher", Path(__file__).with_name("mapped_artifact_publish.py"))
    publisher = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(publisher)
    binding = motif_binding(config_file, source_root, base)
    lineage = json.loads(Path(lineage_file).read_text())
    require(lineage.get("schema") == MOTIF_SCHEMA and lineage.get("scientific_input") == binding
            and lineage.get("provenance"), "Preserved native motif lineage differs from the complete base/producer")
    config, source_root, paths = settings(config_file, source_root, base)
    args = SimpleNamespace(artifact=Path(artifact), mapping=Path(mapping), input=inputs,
                           cache_root=Path(config["cache_root"]), output_dir=Path(artifact).parent,
                           assembly=HEADER.assembly_name(config["assembly"]), config=Path(config_file),
                           producer=source_root / "scripts/am_pick_missense_constrained_segments.py", kind="motifs",
                           annotation_base=paths["base"], native_lineage=Path(lineage_file))
    publisher.adopt_native(args, Path(lineage["artifact_path"]), lineage)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    commands = parser.add_subparsers(dest="command", required=True)
    for name in ("path", "capture", "check", "recover", "publish", "restore", "adopt", "runtime-associated", "motif-binding", "motif-check", "adopt-motif", "motif-runtime"):
        command = commands.add_parser(name)
        command.add_argument("--config", type=Path, required=name != "motif-runtime")
        command.add_argument("--source-root", type=Path, required=name != "motif-runtime")
        command.add_argument("--base", type=Path)
        if name == "publish":
            command.add_argument("--staged-converted", type=Path, required=True)
            command.add_argument("--staged-annotation", type=Path, required=True)
            command.add_argument("--expected-inputs", type=Path, required=True)
        if name == "capture":
            command.add_argument("--output", type=Path, help="Write the captured metadata beside the private scientific stage")
        if name in ("adopt", "adopt-motif"):
            command.add_argument("--lineage", type=Path, required=True)
        if name in ("motif-check", "adopt-motif"):
            command.add_argument("--artifact", type=Path, required=True)
        if name == "adopt-motif":
            command.add_argument("--mapping", type=Path, required=True)
            command.add_argument("--input", type=Path, action="append", default=[])
    args = parser.parse_args()
    common = (args.config, args.source_root, args.base)
    try:
        if args.command == "motif-runtime":
            print(json.dumps(motif_runtime(), sort_keys=True))
        elif args.command == "path":
            print(settings(*common)[2]["base"])
        elif args.command == "capture":
            value = capture(*common)
            if args.output:
                atomic_json(args.output, value)
            else:
                print(json.dumps(value, sort_keys=True))
        elif args.command == "check":
            check(*common)
        elif args.command == "recover":
            recover(*common)
        elif args.command == "publish":
            publish(args.config, args.source_root, args.staged_converted, args.staged_annotation, args.expected_inputs, args.base)
        elif args.command == "restore":
            restore(*common)
        elif args.command == "adopt":
            adopt(args.config, args.source_root, args.lineage, args.base)
        elif args.command == "runtime-associated":
            runtime_associated(*common)
        elif args.command == "motif-binding":
            print(json.dumps(motif_binding(*common), sort_keys=True))
        elif args.command == "adopt-motif":
            adopt_motif(args.config, args.source_root, args.artifact, args.lineage, args.mapping, args.input, args.base)
        else:
            motif_check(args.config, args.source_root, args.artifact, args.base)
    except (OSError, ValueError, KeyError, TypeError, subprocess.SubprocessError) as error:
        parser.exit(1, "AlphaMissense complete annotation: " + str(error) + "\n")


if __name__ == "__main__":
    main()
