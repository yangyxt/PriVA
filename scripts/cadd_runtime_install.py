#!/usr/bin/env python3
"""Install or verify pinned CADD runtimes/resources for annotate, calculate or both.

This route never runs VEP or prediction and never activates a scientific profile.
Calculation uses the pinned original Python through raw scores and PHRED.
Resource presence does not qualify shared PriVA VEP features.
"""
import argparse
from contextlib import contextmanager
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


HERE = Path(__file__).resolve().parent


def sibling(name):
    spec = importlib.util.spec_from_file_location("priva_" + name, HERE / (name + ".py"))
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


runtime = sibling("cadd_runtime")
cache = sibling("cadd_cache_install")
transcripts = sibling("cadd_transcript_resources")
native = sibling("cadd_native")
SOURCE_COMMIT = "1fead0c1e02d3f0a8599bdb04e6b3a2dd70b4866"
PYTHON_CALCULATION_SHA256 = "18d10c4533fc9275a12d0ac675b55903b2ee835d4c6f31de2a3a86452c26382d"
SOURCE_PATHS = ("src", "config", "envs", "Snakefile", "data/models", "CADD.sh")


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def sha256(path):
    return cache.checksum(Path(path), "sha256")


def run(command, **kwargs):
    print("+ " + " ".join(map(str, command)), flush=True)
    return subprocess.run(list(map(str, command)), check=True, **kwargs)


def load_pins(path):
    pins = json.loads(path.read_text())
    require(pins.get("schema_version") == 1 and pins.get("cadd_version") == "v1.7",
            "This installer requires the reviewed v1.7 manifest")
    require(pins.get("source_commit") == SOURCE_COMMIT,
            "The CADD source pin must remain the reviewed v1.7.2 revision")
    require("standalone" not in pins,
            "The Rust CADD branch is suspended; select the original-Python runtime manifest")
    require(pins.get("calculation", {}).get("implementation") == "all_four_calculation_stages_upstream_python"
            and pins["calculation"].get("wrapper_sha256") == PYTHON_CALCULATION_SHA256
            and sha256(HERE / "cadd_python_calculate.sh") == PYTHON_CALCULATION_SHA256,
            "The original Python calculation command sequence must retain its reviewed pin")
    require(re.fullmatch(r"docker://[^\s]+@sha256:[0-9a-f]{64}", pins.get("image", "")),
            "CADD image must use an immutable OCI digest")
    require(isinstance(pins.get("source_url"), str) and pins["source_url"].startswith("https://"),
            "Source URL must use HTTPS")
    return pins


def checkout(path, url, revision, paths, verify_only=False):
    """Reuse exactly pinned tracked bytes; never clean or reset an existing checkout."""
    if not path.exists():
        require(not verify_only, f"CADD source is not installed: {path}")
        path.parent.mkdir(parents=True, exist_ok=True)
        run(["git", "clone", "--no-checkout", url, path])
        run(["git", "-C", path, "checkout", "--detach", revision])
    actual = subprocess.check_output(["git", "-C", str(path), "rev-parse", "HEAD"], text=True).strip()
    require(actual == revision, f"Pinned source revision differs: {path}; existing checkout retained")
    run(["git", "-C", path, "diff", "--quiet", "HEAD", "--", *paths])
    untracked = subprocess.check_output(
        ["git", "-C", str(path), "ls-files", "--others", "--exclude-standard", "--", *paths],
        text=True).strip()
    require(not untracked, f"Untracked files can shadow pinned source in {path}; existing files retained")
    return path


@contextmanager
def locked(path):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("a") as lock:
        try:
            fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        except BlockingIOError as error:
            raise RuntimeError(f"CADD installation writer is active: {path}; retry after it finishes") from error
        yield


def publish_json(path, record):
    """Publish a new immutable receipt, or reuse an identical prior one."""
    path.parent.mkdir(parents=True, exist_ok=True)
    if os.path.lexists(path):
        previous = json.loads(path.read_text())
        for key, value in record.items():
            if key not in ("checked_utc", "slurm_job_id"):
                require(previous.get(key) == value, f"Existing CADD receipt differs: {path}; retained")
        return path
    temporary = None
    try:
        with tempfile.NamedTemporaryFile(mode="w", prefix="." + path.name + ".", dir=path.parent,
                                         delete=False) as output:
            json.dump(record, output, indent=2, sort_keys=True)
            output.write("\n")
            temporary = Path(output.name)
        os.link(temporary, path)
    finally:
        if temporary is not None:
            temporary.unlink(missing_ok=True)
    return path


def image_and_routes(folder, pins, state, temporary, verify_only=False):
    folder.mkdir(parents=True, exist_ok=True)
    image = folder / "cadd-upstream-0.1.1.sif"
    origin = folder / "image.origin.json"
    environment = dict(os.environ)
    environment["APPTAINER_CACHEDIR"] = str(folder / "image-cache")
    environment["APPTAINER_TMPDIR"] = str(temporary / "cadd-apptainer")
    environment["APPTAINERENV_MPLCONFIGDIR"] = str(temporary / "matplotlib")
    environment["OMP_NUM_THREADS"] = "1"
    if not image.exists():
        require(not verify_only, f"CADD runtime image is not installed: {image}")
        require(not origin.exists(), "CADD image is absent but its origin receipt exists; recover explicitly")
        Path(environment["APPTAINER_TMPDIR"]).mkdir(parents=True, exist_ok=True)
        partial = image.with_name(image.name + ".partial")
        # A SIF conversion is not a byte-range transfer. Retain interrupted SIFs
        # and require explicit recovery rather than silently replacing them.
        require(not os.path.lexists(partial), f"Retained incomplete CADD image: {partial}; recover explicitly")
        run(["apptainer", "pull", partial, pins["image"]], env=environment)
        digest = sha256(partial)
        os.link(partial, image)
        partial.unlink()
        publish_json(origin, {"image_source": pins["image"], "image_sha256": digest})
    else:
        previous_path = origin if origin.exists() else state / "cadd-runtime.ready.json"
        require(previous_path.is_file(), "Existing CADD image has no verified source receipt")
        previous = json.loads(previous_path.read_text())
        digest = sha256(image)
        require(previous.get("image_source") == pins["image"] and previous.get("image_sha256") == digest,
                "Existing CADD image differs from its verified source receipt; bytes retained")
    before = cache.identity(image)
    if verify_only:
        wrapper = folder / "bin/cadd-container-exec"
        require(runtime.container_wrapper_matches(wrapper), "CADD runtime wrapper is absent or differs")
        for name in runtime.CADD_ENTRY_POINTS:
            alias = folder / "bin" / name
            require(alias.is_symlink() and os.readlink(alias) == wrapper.name,
                    f"CADD named interpreter route is absent or differs: {alias}")
        routes = {"path": str(wrapper), "sha256": sha256(wrapper),
                  "entry_points": {name: str(folder / "bin" / name) for name in runtime.CADD_ENTRY_POINTS}}
    else:
        routes = runtime.install_cadd_entry_points(folder)
    return {"path": str(image), "sha256": digest, "identity": before,
            "source": pins["image"], "wrappers": routes}, environment


def annotation_resources(config, assembly, source, cache_root, verify_only=False, receipt_root=None):
    version = config.get("cadd_version")
    require(version == "v1.7", "CADD resources must match the reviewed v1.7 model")
    prefix = f"cadd_{assembly}_anno_"
    url, md5 = config.get(prefix + "url"), config.get(prefix + "md5")
    require(isinstance(url, str) and re.match(r"^(?:https?|file)://", url), f"Set {prefix}url")
    require(isinstance(md5, str) and re.fullmatch(r"[0-9a-fA-F]{32}", md5), f"Set {prefix}md5")
    resources = {"assembly": assembly, "version": version, "annotation": (url, md5.lower())}
    folder = source / "data/annotations"
    target = folder / f"{assembly}_{version}"
    archive = cache.annotation_archive(resources, folder, cache_root)
    member_receipt = Path(receipt_root) / (assembly + ".members.ready.json") if receipt_root else target / cache.RECEIPT
    # An explicitly selected legacy source can retain its already verified
    # archive while new runtime/package caches use a different configured root.
    # Never bypass an existing canonical archive or acquire bytes outside it.
    if receipt_root and not os.path.lexists(archive):
        retained = source.parent.parent / ".downloads" / archive.name
        if retained.is_file():
            archive = retained
    if verify_only or receipt_root is not None:
        require(member_receipt.is_file(),
                f"CADD annotation member verification is incomplete: {member_receipt}; run installation first")
        previous = json.loads(member_receipt.read_text())
        require(all(previous.get(key) == expected for key, expected in
                    (("schema", 1), ("assembly", assembly), ("version", version),
                     ("archive_url", url), ("archive_md5", md5.lower()))),
                "CADD annotation receipt does not match configured assembly/source/checksum")
        require(archive.is_file(), f"Verified CADD archive is missing: {archive}")
        require(cache.identity(archive) == previous.get("archive_identity")
                or cache.checksum(archive) == md5.lower(), "Existing CADD annotation archive changed")
        cache.verify_receipt(target, previous)
    else:
        # This is the same lock held by frozen GRCh38 deployment 1515113.
        with locked(source / "data/.priva-cadd.install.lock"):
            cache.annotations(resources, folder, archive_path=archive)
        previous = json.loads(member_receipt.read_text())
    return {"path": str(target), "archive": str(archive), "archive_md5": md5.lower(),
            "archive_identity": cache.identity(archive), "member_receipt": str(member_receipt),
            "member_receipt_sha256": sha256(member_receipt),
            "verified_members": len(previous["members"]),
            "integrity": "publisher_md5_and_extracted_member_sha256"}


def calculation_resources(source, assembly):
    names = [f"config/impute_{assembly}_v1.7.cfg",
             f"data/models/{assembly}/CADDv1.7-{assembly}.mod",
             f"data/models/{assembly}/conversionTable_CADDv1.7-{assembly}.txt",
             *(f"src/scripts/{name}.py" for name in
               ("trackTransformation", "predictSKmodel", "max_line_hierarchy", "appendPHREDscore"))]
    result = []
    for name in names:
        path = source / name
        require(path.is_file() and path.stat().st_size > 0, f"Missing original calculation asset: {path}")
        result.append({"path": str(path), "sha256": sha256(path), "identity": cache.identity(path)})
    return result


def python_calculation(program, source, folder, assembly, verify_only=False):
    route = runtime.install_python_calculation(program, verify_only)
    route["environment"] = {
        "CADD_SOURCE_ROOT": str(source), "CADD_PYTHON2": str(folder / "bin/cadd-python2"),
        "CADD_IMPUTE_CONFIG": str(source / f"config/impute_{assembly}_v1.7.cfg"),
        "CADD_MODEL": str(source / f"data/models/{assembly}/CADDv1.7-{assembly}.mod"),
        "CADD_PHRED_TABLE": str(source / f"data/models/{assembly}/conversionTable_CADDv1.7-{assembly}.txt"),
    }
    route["upstream_source_commit"] = SOURCE_COMMIT
    return route


def recognized_legacy_runtime(folder, state, pins):
    """Select an existing legacy runtime only with its explicit immutable origin."""
    image = folder / "cadd-upstream-0.1.1.sif"
    if not image.is_file():
        return False
    for receipt in (folder / "image.origin.json", state / "cadd-runtime.ready.json"):
        if receipt.is_file():
            record = json.loads(receipt.read_text())
            digest = record.get("image_sha256")
            if (record.get("image_source") == pins["image"] and isinstance(digest, str)
                    and re.fullmatch(r"[0-9a-f]{64}", digest)):
                return True
    raise RuntimeError(f"Existing legacy CADD image lacks a matching origin receipt: {folder}; "
                       "set --runtime-dir and verified --state-dir explicitly; existing bytes retained")


def installation_paths(config, args, pins=None):
    """Anchor legacy configurations to their selected source, never an empty root."""
    configured_root = config.get("cache_root")
    require(configured_root is None or isinstance(configured_root, str), "cache_root must be a path string")
    cache_root = Path(configured_root.strip()).absolute() if configured_root and configured_root.strip() else None
    selected_source = args.source or config.get("cadd_base_dir")
    if isinstance(selected_source, str):
        selected_source = selected_source.strip() or None
    if selected_source is None:
        require(cache_root is not None, "Set cache_root, cadd_base_dir or --source before CADD runtime installation")
        selected_source = cache_root / "CADD/CADD-scripts-1.7.2"
    require(isinstance(selected_source, (str, Path)), "cadd_base_dir must be a path string")
    source = Path(selected_source).absolute()
    folder = args.runtime_dir
    if folder is None:
        folder = cache_root / "runtime/cadd-1.7.2" if cache_root else source / ".priva-runtime"
    state = args.state_dir or folder / "receipts"
    # Explicit CLI/environment scratch wins; main_install's configured tmp_dir
    # precedes the historical fallback for configurations without that key.
    temporary = args.tmp_root or config.get("tmp_dir") or state / "tmp"
    program = args.calculation_program or folder / "bin/cadd-python-calculate"
    return (cache_root, *(Path(path).absolute() for path in
                          (source, folder, state, temporary, program)))


def select_backend(args, folder):
    """Honor an explicit choice; only reuse a container explicitly selected by path."""
    if args.backend != "auto":
        return args.backend
    if args.runtime_dir is not None and (folder / "cadd-upstream-0.1.1.sif").is_file():
        return "container"
    return "native"


def native_environment_plan(folder, source, name, pins):
    """Describe original YAML bytes and measured scientific package constraints."""
    definition = pins["native"]["environments"][name]
    env_file = source / definition["environment_file"]
    plan = {"schema_version": 1, "name": name, "prefix": str(folder / "envs" / name),
            "source_commit": pins["source_commit"], "environment_file": str(env_file),
            "environment_sha256": sha256(env_file), "specs": definition.get("specs", []),
            "channels": definition.get("channels", [])}
    if definition.get("constraints"):
        plan["constraints"] = definition["constraints"]
    return plan


def native_package_manager(selected=None):
    """Keep the bootstrap's initialized wrapper, then discover an existing manager."""
    if selected:
        return str(selected)
    candidates = (os.environ.get("MAMBA_EXE"), shutil.which("micromamba"),
                  shutil.which("mamba"), os.environ.get("CONDA_EXE"), shutil.which("conda"))
    for candidate in candidates:
        if (candidate and Path(candidate).name.startswith(("micromamba", "mamba", "conda"))
                and os.access(candidate, os.X_OK) and Path(candidate).is_file()):
            return str(candidate)
    raise RuntimeError("Native CADD installation requires the existing micromamba, mamba or conda executable")


def native_install_commands(manager, operation, folder, prefix, definition, staged):
    """Use each manager's supported CLI, preserving upstream Conda/pip ownership."""
    kind = Path(manager).name
    if kind.startswith("micromamba"):
        command = [manager, operation, "--yes", "--no-rc", "--no-env", "--root-prefix",
                   str(folder / "mamba"), "--prefix", str(prefix)]
        if definition.get("specs"):
            command += ["--override-channels", "--channel-priority", "flexible"]
            for channel in definition["channels"]:
                command += ["--channel", channel]
            command += definition["specs"]
        else:
            command += ["--file", str(staged), *definition.get("constraints", [])]
        return command, None
    require(kind.startswith(("conda", "mamba")), f"Unrecognized native CADD package manager: {manager}")
    command = [manager, operation, "--yes", "--prefix", str(prefix), "--override-channels"]
    if operation == "create":
        command += ["--no-default-packages"]
    pip_specs = []
    if definition.get("specs"):
        channels, specs = definition["channels"], definition["specs"]
    else:
        # Conda versions predating environment-specifier plugins cannot accept
        # YAML and additional CLI constraints with create/install. Read only the
        # unchanged upstream channel/dependency lists, then run its pip phase in
        # the selected prefix, as conda env create itself does.
        import yaml
        document = yaml.safe_load(staged.read_text())
        channels, specs = document["channels"], []
        for dependency in document["dependencies"]:
            if isinstance(dependency, str):
                specs.append(dependency)
            else:
                require(isinstance(dependency, dict) and set(dependency) == {"pip"}
                        and isinstance(dependency["pip"], list), "Unsupported original CADD dependency section")
                pip_specs += dependency["pip"]
        specs += definition.get("constraints", [])
    for channel in channels:
        if channel != "nodefaults":
            command += ["--channel", channel]
    command += specs
    pip_command = None
    if pip_specs:
        requirements = staged.with_suffix(".pip.txt")
        payload = "\n".join(pip_specs) + "\n"
        if requirements.exists():
            require(not requirements.is_symlink() and requirements.read_text() == payload,
                    f"Existing native CADD pip requirements differ: {requirements}; retained")
        else:
            with requirements.open("x") as stream:
                stream.write(payload)
        pip_command = [str(prefix / "bin/python"), "-m", "pip", "install", "--requirement", str(requirements)]
    return command, pip_command


def native_environments(folder, source, names, temporary, pins, verify_only=False, manager=None):
    """Install original Python environments without Snakemake or VEP execution.

    The source-owned YAML specifies feature environments. Core package versions
    retain the measured reference. Completed prefixes are never solved again;
    incomplete installer-owned prefixes resume using install, retaining bytes.
    """
    environment = dict(os.environ)
    environment.pop("PYTHONHOME", None)
    environment.pop("PYTHONPATH", None)
    environment.update(MAMBA_ROOT_PREFIX=str(folder / "mamba"),
                       CONDA_PKGS_DIRS=str(folder / "mamba/pkgs"),
                       CONDA_ENVS_PATH=str(folder / "envs"), CONDA_CHANNEL_PRIORITY="flexible",
                       PIP_CACHE_DIR=str(folder / "pip-cache"),
                       XDG_CACHE_HOME=str(folder / "runtime-cache"),
                       MPLCONFIGDIR=str(folder / "runtime-cache/matplotlib"),
                       KERAS_HOME=str(folder / "runtime-cache/keras"),
                       TORCH_HOME=str(folder / "runtime-cache/torch"),
                       KIPOI_ENV_DB_PATH=str(folder / "runtime-cache/kipoi-envs.json"),
                       TMPDIR=str(temporary), PYTHONDONTWRITEBYTECODE="1",
                       PYTHONNOUSERSITE="1", OMP_NUM_THREADS="1",
                       OPENBLAS_NUM_THREADS="1", MKL_NUM_THREADS="1", NUMEXPR_NUM_THREADS="1")
    state = folder / "native-state"
    state.mkdir(parents=True, exist_ok=True)
    records = {}
    for name in names:
        definition = pins["native"]["environments"][name]
        prefix = folder / "envs" / name
        env_file = source / definition["environment_file"]
        plan = native_environment_plan(folder, source, name, pins)
        planned = state / (name + ".plan.json")
        ready = state / (name + ".ready.json")
        if ready.is_file():
            record = json.loads(ready.read_text())
            # Authorized cache relocation changes this descriptive source path,
            # while commit, exact YAML bytes and package constraints stay fixed.
            identity = lambda value: {key: item for key, item in value.items() if key != "environment_file"}
            require(identity(record.get("plan", {})) == identity(plan),
                    f"Native CADD environment plan differs: {prefix}; retained")
            runtime.verify_native_environment(prefix, record)
        else:
            require(not verify_only, f"Native CADD environment is not installed: {prefix}")
            if prefix.exists() and not planned.exists():
                raise RuntimeError(f"Unmanaged existing CADD environment retained: {prefix}")
            publish_json(planned, plan)
            manager = native_package_manager(manager)
            # create may delete an existing prefix. Only use it for an absent
            # prefix; an interrupted owned prefix is resumed with install.
            if prefix.exists() and not (prefix / "conda-meta/history").is_file():
                retained = prefix.with_name(prefix.name + ".incomplete." + datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%S%f"))
                prefix.rename(retained)
                print(f"Retained incomplete native CADD prefix at {retained}", flush=True)
            operation = "install" if prefix.exists() else "create"
            staged = None
            if not definition.get("specs"):
                # Micromamba writes transient pip requirements beside the YAML.
                # Stage identical bytes inside the runtime so source stays clean.
                staged = state / "environment-files" / (name + ".yml")
                staged.parent.mkdir(exist_ok=True)
                if staged.exists():
                    require(not staged.is_symlink() and staged.read_bytes() == env_file.read_bytes(),
                            f"Existing staged CADD environment differs: {staged}; retained")
                else:
                    with staged.open("xb") as stream:
                        stream.write(env_file.read_bytes())
            command, pip_command = native_install_commands(manager, operation, folder, prefix, definition, staged)
            run(command, env=environment)
            if pip_command:
                pip_environment = dict(environment, CONDA_PREFIX=str(prefix),
                                       PATH=str(prefix / "bin") + os.pathsep + environment.get("PATH", ""))
                run(pip_command, env=pip_environment)
            require((prefix / "bin/python").is_file() and (prefix / "conda-meta/history").is_file(),
                    f"Native CADD package manager did not complete {prefix}")
            metadata = {p.name: sha256(p) for p in sorted((prefix / "conda-meta").glob("*.json"))}
            record = {"plan": plan, "conda_metadata_sha256": metadata,
                      "package_manager": str(manager), "command": list(map(str, command)),
                      "status": "installed_not_scientifically_qualified"}
            if pip_command:
                record["pip_command"] = pip_command
            publish_json(ready, record)
        records[name] = {"prefix": str(prefix), "receipt": str(ready), "sha256": sha256(ready)}
    routes = runtime.install_cadd_entry_points(folder, native=True, names=names)
    return {"backend": "native", "wrappers": routes, "environments": records}, environment


def install(config, args, pins):
    assembly = {"hg19": "GRCh37", "GRCh37": "GRCh37", "hg38": "GRCh38", "GRCh38": "GRCh38"}.get(config.get("assembly"))
    require(assembly is not None, "assembly must be hg19/GRCh37 or hg38/GRCh38")
    stages = ("annotate", "calculate") if args.stage == "both" else (args.stage,)
    cache_root, source, folder, state, temporary, program = installation_paths(config, args, pins)
    backend = select_backend(args, folder)
    if backend == "container":
        require(shutil.which("apptainer"), "Install Apptainer before selecting the CADD container runtime")
    require(shutil.which("git"), "Install Git before CADD runtime installation")
    state.mkdir(parents=True, exist_ok=True)
    temporary.mkdir(parents=True, exist_ok=True)
    with locked(state / "cadd-runtime.install.lock"):
        checkout(source, pins["source_url"], pins["source_commit"], SOURCE_PATHS, args.verify_only)
        names = ("cadd-python2",)
        if "annotate" in stages:
            names += ("cadd-esm-python", "cadd-regseq-python")
            if assembly == "GRCh38":
                names += ("cadd-mmsplice-python",)
        if backend == "native":
            image, environment = native_environments(folder, source, names, temporary, pins,
                                                     args.verify_only, args.package_manager)
        else:
            image, environment = image_and_routes(folder, pins, state, temporary, args.verify_only)
        annotation = None
        selected_transcripts = None
        if "annotate" in stages:
            annotation = annotation_resources(config, assembly, source, cache_root, args.verify_only,
                                              args.annotation_member_state)
            selected_transcripts = transcripts.install(config, verify_only=args.verify_only)
        check_options = {"names": names, "include_vep": False}
        if backend == "native":
            check_options["native"] = True
        environments = runtime.check_cadd_environments(folder, environment, **check_options)
        assets = runtime.cadd_asset_readiness(source, args.asset_view_dir or folder, assemblies=(assembly,), stages=stages,
                                              hub_assembly=assembly)
        if "annotate" in stages:
            require(assets["annotation_assets"][assembly]["status"] != "incomplete",
                    f"Required {assembly} feature-model assets are incomplete")
        calculation = calculation_resources(source, assembly) if "calculate" in stages else []
        calculation_runtime = python_calculation(program, source, folder, assembly, args.verify_only) if "calculate" in stages else None
        if backend == "container":
            require(cache.identity(Path(image["path"])) == image["identity"], "CADD runtime image changed during verification")
        checkout(source, pins["source_url"], pins["source_commit"], SOURCE_PATHS, True)
        record = {"schema_version": 1, "stage": "cadd-install", "selection": args.stage,
                  "status": "installed_not_scientifically_qualified", "assembly": assembly,
                  "cadd_version": pins["cadd_version"], "source": str(source),
                  "source_commit": pins["source_commit"], "backend": backend, "runtime": image,
                  "image": image if backend == "container" else None,
                  "environments": environments, "assets": assets, "annotation": annotation,
                  "selected_transcript_candidate": selected_transcripts,
                  "calculation_assets": calculation, "calculation_runtime": calculation_runtime,
                  "manifest": {"path": str(args.manifest.absolute()), "sha256": sha256(args.manifest)},
                  "producer": {"path": str(Path(__file__).resolve()), "sha256": sha256(Path(__file__)),
                               "runtime_helper_sha256": sha256(HERE / "cadd_runtime.py"),
                               "cache_helper_sha256": sha256(HERE / "cadd_cache_install.py")},
                  "qualification": {"resource_runtime_installation": True, "vep_executed": False,
                                    "model_inference_performed": False, "shared_vep_qualified": False,
                                    "scientific_profile_activated": False,
                                    "legacy_transcript_release": 110,
                                    "legacy_transcript_release_evidence": {
                                        "config_path": str(source / f"config/config_{assembly}_v1.7.yml"),
                                        "config_key": "EnsemblDB",
                                        "esm_peptide_path": str(source / f"data/annotations/{assembly}_v1.7/esm/pep.110.fa")},
                                    "legacy_coordinate_annotation_release": 75 if assembly == "GRCh37" else 110,
                                    "legacy_coordinate_annotation_evidence": {
                                        "config_path": str(source / f"config/references_{assembly}_v1.7.cfg"),
                                        "config_keys": ["Exon/Dst2Splice", "Transcript/MinDistTSX"],
                                        "exon_asset_path": str(source / ("data/annotations/GRCh37_v1.7/reference/annotation_exons_v75.gtf.gz" if assembly == "GRCh37" else "data/annotations/GRCh38_v1.7/exons/annotation_exons_GRCh38.110.gtf.gz")),
                                        "transcript_asset_path": str(source / ("data/annotations/GRCh37_v1.7/reference/annotation_transcripts_v75.tsv.gz" if assembly == "GRCh37" else "data/annotations/GRCh38_v1.7/transcripts/annotation_transcripts_GRCh38.110.tsv.gz"))},
                                    "shared_priva_vep_release": str(config.get("vep_installed_version", ""))}}
        fingerprint = hashlib.sha256(json.dumps(record, sort_keys=True).encode()).hexdigest()
        record.update({"fingerprint": fingerprint, "checked_utc": datetime.now(timezone.utc).isoformat(),
                       "slurm_job_id": os.environ.get("SLURM_JOB_ID")})
        receipt = state / f"cadd-install.{assembly}.{args.stage}.{fingerprint[:16]}.ready.json"
        publish_json(receipt, record)
        deployment = (native.install_profile(config, receipt, args.verify_only)
                      if backend == "native" and args.stage == "both" else None)
        print(json.dumps({"status": record["status"], "assembly": assembly, "selection": args.stage,
                          "receipt": str(receipt), "native_profile": deployment,
                          "scientific_profile_activated": False}), flush=True)
        return receipt


def arguments(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    config = parser.add_mutually_exclusive_group(required=True)
    config.add_argument("--config", type=Path)
    config.add_argument("--config-json")
    parser.add_argument("--stage", choices=("annotate", "calculate", "both"), default="both")
    parser.add_argument("--manifest", type=Path, default=HERE / "cadd_runtime_pins.json")
    parser.add_argument("--source-only", action="store_true",
                        help="Prepare the reviewed pinned Git source for CADD_install; acquire no data/image")
    parser.add_argument("--backend", choices=("auto", "native", "container"), default="auto",
                        help="Fresh installs use native Python; an explicitly selected existing image retains its container route")
    parser.add_argument("--package-manager", help="Existing native micromamba, mamba or conda executable (defaults to initialized manager/PATH)")
    for flag, variable in (("runtime-dir", "PRIVA_CADD_RUNTIME_DIR"), ("state-dir", "PRIVA_CADD_STATE_DIR"),
                           ("tmp-root", "PRIVA_CADD_TMP_ROOT"), ("source", "PRIVA_CADD_SOURCE"),
                           ("calculation-program", "PRIVA_CADD_CALCULATION_PROGRAM")):
        value = os.environ.get(variable)
        parser.add_argument("--" + flag, type=Path, default=Path(value) if value else None)
    parser.add_argument("--asset-view-dir", type=Path, help="Writable runtime view for lightweight ESM checkpoint links")
    parser.add_argument("--annotation-member-state", type=Path,
                        help="Reuse completed external member receipts without writing the selected legacy annotation tree")
    parser.add_argument("--verify-only", action="store_true", help="Require installed assets; do not acquire missing source/image/resources")
    return parser.parse_args(argv)


def main(argv=None):
    args = arguments(argv)
    try:
        if args.config_json is not None:
            config = json.loads(args.config_json)
        else:
            import yaml
            config = yaml.safe_load(args.config.read_text())
        require(isinstance(config, dict), "CADD configuration must be a mapping")
        pins = load_pins(args.manifest)
        if args.source_only:
            require(args.source is not None, "--source-only requires --source")
            require(config.get("cadd_zip_download_url") == pins.get("source_archive_url"),
                    "Configured CADD archive is not the reviewed v1.7.2 source; "
                    "custom archive/source bytes are retained for explicit review")
            checkout(args.source, pins["source_url"], pins["source_commit"], SOURCE_PATHS, args.verify_only)
            print(json.dumps({"source": str(args.source), "source_commit": pins["source_commit"],
                              "status": "pinned_source_prepared"}), flush=True)
        else:
            install(config, args, pins)
    except (OSError, RuntimeError, ValueError, KeyError, subprocess.SubprocessError) as error:
        print(f"ERROR: CADD runtime installation failed: {error}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
