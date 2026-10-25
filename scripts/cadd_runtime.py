#!/usr/bin/env python3
"""Shared CADD image routes and bounded checks, preserving the deployed v1.7.2 runtime.

Installation integrity is separate from shared-VEP scientific qualification.
This module is also re-exported by deploy_cluster_resources.py for compatibility.
"""
import hashlib
import json
import os
from pathlib import Path
import re
import subprocess
import tempfile

def checksum(path, algorithm="sha256"):
    digest = hashlib.new(algorithm)
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def verify_native_environment(prefix, record):
    """Read-only reuse guard shared by installation and native execution."""
    prefix = Path(prefix)
    if not (prefix / "bin/python").is_file() or not (prefix / "conda-meta/history").is_file():
        raise RuntimeError(f"Completed native CADD environment is incomplete: {prefix}")
    for member, expected in record["conda_metadata_sha256"].items():
        path = prefix / "conda-meta" / member
        if not path.is_file():
            raise RuntimeError(f"Native CADD package metadata changed: {path}; retained")
        before = path.stat()
        measured = checksum(path)
        after = path.stat()
        fields = ("st_dev", "st_ino", "st_size", "st_mtime_ns", "st_ctime_ns")
        if measured != expected or any(getattr(before, key) != getattr(after, key) for key in fields):
            raise RuntimeError(f"Native CADD package metadata changed: {path}; retained")
    if set(record["conda_metadata_sha256"]) != {p.name for p in (prefix / "conda-meta").glob("*.json")}:
        raise RuntimeError(f"Native CADD package membership changed: {prefix}; retained")



CADD_ENVIRONMENTS = {
    "cadd-python2": ("/conda-envs/fe7049a03d989a2faa6d7688f5987d09", 2,
                     ("sklearn", "numpy", "pandas", "scipy", "pysam")),
    "cadd-esm-python": ("/conda-envs/10a710eb7e199e21064fc58a16f9ea4f", 3,
                        ("torch", "esm", "Bio", "click")),
    "cadd-regseq-python": ("/conda-envs/f675d936396bfaffdd96e0282b787268", 3,
                           ("tensorflow", "numpy", "seqiolib", "vcfpy", "pyfaidx", "pybedtools")),
    "cadd-mmsplice-python": ("/conda-envs/27e6f085c62573d05046d5330f07505d", 3,
                            ("tensorflow", "keras", "mmsplice", "concise", "numpy", "pandas", "tqdm")),
}
CADD_ENTRY_POINTS = (*CADD_ENVIRONMENTS, "cadd-vep")
ESM_MODELS = tuple(f"esm1v_t33_650M_UR90S_{number}" for number in range(1, 6))
LEGACY_CONTAINER_WRAPPER_SHA256 = "b28824617dbb6c0df12a441162c5abfc2d93b63d78746b92f729f020ebd50827"


def container_wrapper_matches(path):
    """Recognize the portable wrapper or the retained, measured HGSC reference."""
    return (path.is_file() and not path.is_symlink() and os.access(path, os.X_OK)
            and (path.read_bytes() == Path(__file__).with_name("cadd_container_exec.sh").read_bytes()
                 or checksum(path) == LEGACY_CONTAINER_WRAPPER_SHA256))


def install_cadd_entry_points(folder, native=False, names=None):
    """Publish only the reviewed wrapper and relative aliases; reuse is a no-op."""
    wrapper_source = Path(__file__).with_name("cadd_native_exec.sh" if native else "cadd_container_exec.sh")
    payload = wrapper_source.read_bytes()
    binary = folder / "bin"
    if binary.is_symlink():
        raise RuntimeError(f"Refusing a symlinked CADD runtime bin directory: {binary}")
    binary.mkdir(exist_ok=True)
    wrapper = binary / ("cadd-native-exec" if native else "cadd-container-exec")
    selected = tuple(names) if names is not None else CADD_ENTRY_POINTS
    if (not native and wrapper.exists() and container_wrapper_matches(wrapper)):
        payload = wrapper.read_bytes()  # Keep an explicitly selected reference immutable.
    if wrapper.is_symlink() or (wrapper.exists() and wrapper.read_bytes() != payload):
        raise RuntimeError(f"Existing CADD wrapper differs from the reviewed source: {wrapper}")
    for name in selected:
        target = binary / name
        if target.is_symlink():
            if os.readlink(target) != wrapper.name:
                raise RuntimeError(f"Existing CADD entry point has a different route: {target}")
        elif target.exists():
            raise RuntimeError(f"Existing CADD entry point is not the expected symlink: {target}")
    if not wrapper.exists():
        with wrapper.open("xb") as stream:
            stream.write(payload)
        wrapper.chmod(0o755)
    if not os.access(wrapper, os.X_OK):
        raise RuntimeError(f"CADD wrapper is not executable: {wrapper}")
    for name in selected:
        target = binary / name
        if not target.is_symlink():
            target.symlink_to(wrapper.name)
    return {"path": str(wrapper), "sha256": checksum(wrapper),
            "entry_points": {name: str(binary / name) for name in selected}}


def install_python_calculation(program, verify_only=False):
    """Publish the unchanged upstream-Python command route, retaining other files."""
    source = Path(__file__).with_name("cadd_python_calculate.sh")
    payload = source.read_bytes()
    program = Path(program)
    if program.is_symlink():
        raise RuntimeError(f"CADD Python calculation program must be a regular frozen file: {program}")
    if program.exists():
        if not program.is_file() or program.read_bytes() != payload or not os.access(program, os.X_OK):
            raise RuntimeError(f"Existing CADD Python calculation route differs: {program}; retained")
    elif verify_only:
        raise RuntimeError(f"CADD Python calculation route is not installed: {program}")
    else:
        program.parent.mkdir(parents=True, exist_ok=True)
        temporary = None
        try:
            with tempfile.NamedTemporaryFile(prefix="." + program.name + ".", dir=program.parent,
                                             delete=False) as handle:
                handle.write(payload)
                temporary = Path(handle.name)
            temporary.chmod(0o755)
            os.link(temporary, program)
        finally:
            if temporary is not None:
                temporary.unlink(missing_ok=True)
    return {"program": str(program), "sha256": checksum(program),
            "implementation": "all_four_calculation_stages_upstream_python",
            "arguments": ["INPUT.annotations.tsv[.gz]", "OUTPUT.scored.tsv", "WORKDIR"]}


def readable_asset(path):
    """Bounded readability, not a substitute for full publisher checks or inference."""
    if not path.is_file() or not os.access(path, os.R_OK) or path.stat().st_size == 0:
        raise RuntimeError(f"Missing, empty or unreadable CADD asset: {path}")
    before = path.stat()
    with path.open("rb") as stream:
        first = stream.read(64)
        stream.seek(max(0, before.st_size - 64))
        last = stream.read(64)
    after = path.stat()
    identity = lambda value: (value.st_dev, value.st_ino, value.st_size, value.st_mtime_ns)
    if identity(before) != identity(after):
        raise RuntimeError(f"CADD asset changed during readability check: {path}")
    return {"path": str(path), "bytes": before.st_size, "mtime_ns": before.st_mtime_ns,
            "sample_sha256": hashlib.sha256(first + last).hexdigest(),
            "check": "first_and_last_64_bytes_readable"}


def cadd_asset_readiness(location, folder, assemblies=("GRCh37", "GRCh38"),
                         stages=("annotate", "calculate"), hub_assembly="GRCh37"):
    models = [readable_asset(location / f"data/models/{assembly}/{name}")
              for assembly in assemblies if "calculate" in stages
              for name in (f"CADDv1.7-{assembly}.mod", f"conversionTable_CADDv1.7-{assembly}.txt")]
    annotations = {}
    for assembly in (assemblies if "annotate" in stages else ()):
        root = location / f"data/annotations/{assembly}_v1.7"
        names = [*(f"esm/{name}.pt" for name in ESM_MODELS), "esm/pep.110.fa",
                 "regseq/Hyperopt400InclNegatives.json", "regseq/Hyperopt400InclNegatives.h5",
                 "regseq/reference.fa", "regseq/reference.fa.fai", "regseq/reference.fa.genome"]
        if assembly == "GRCh38":
            names += ["mmsplice/homo_sapiens.110.gtf", "reference/reference.fa", "reference/reference.fa.fai"]
        present = [readable_asset(root / name) for name in names if (root / name).is_file()]
        missing = [str(root / name) for name in names if not (root / name).is_file()]
        annotations[assembly] = {"assets": present, "missing": missing,
                                 "status": "readable_not_scientifically_qualified" if not missing else "incomplete"}
    # Upstream esmScore sets torch.hub's directory and then resolves named
    # models under checkpoints/. Keep original flat downloads unchanged.
    esm_root = location / f"data/annotations/{hub_assembly}_v1.7/esm"
    hub = folder / ("esm-hub" if hub_assembly == "GRCh37" else f"esm-hub-{hub_assembly}")
    hub_ready = "annotate" in stages and all((esm_root / f"{name}.pt").is_file() for name in ESM_MODELS)
    if hub_ready:
        if hub.is_symlink() or (hub / "checkpoints").is_symlink():
            raise RuntimeError("Refusing a symlinked ESM hub/checkpoints directory")
        checkpoints = hub / "checkpoints"
        checkpoints.mkdir(parents=True, exist_ok=True)
        retargets = []
        for name in ESM_MODELS:
            target, link = esm_root / f"{name}.pt", checkpoints / f"{name}.pt"
            if link.is_symlink():
                previous = os.readlink(link)
                if previous != str(target):
                    try:
                        same_model = link.samefile(target)
                    except OSError:
                        same_model = False
                    if not same_model:
                        raise RuntimeError(f"Existing ESM checkpoint has a different target: {link}")
                    retargets.append((link, target, previous))
            elif link.exists():
                raise RuntimeError(f"Existing ESM checkpoint is not the expected symlink: {link}")
        # A managed-cache handoff can give the same hardlinked model a new
        # pathname. Validate every checkpoint first, then replace only links
        # whose old and selected targets are proven to be the same inode.
        for link, target, previous in retargets:
            with tempfile.TemporaryDirectory(prefix=".cadd-checkpoint-", dir=checkpoints) as temporary:
                updated = Path(temporary) / "checkpoint"
                updated.symlink_to(str(target))
                if (not link.is_symlink() or os.readlink(link) != previous
                        or not link.samefile(target)):
                    raise RuntimeError(f"ESM checkpoint changed during cache handoff: {link}")
                os.replace(updated, link)
            print(f"Retargeted ESM checkpoint to the same original model: {link}: {previous} -> {target}", flush=True)
        for name in ESM_MODELS:
            target, link = esm_root / f"{name}.pt", checkpoints / f"{name}.pt"
            if not link.is_symlink():
                link.symlink_to(target)
    return {"calculation_models": models, "annotation_assets": annotations,
            "esm_hub": {"path": str(hub), "status": "ready" if hub_ready else "unavailable",
                        "checkpoints": [{"link": str(hub / "checkpoints" / f"{name}.pt"),
                                         "target": str(esm_root / f"{name}.pt")} for name in ESM_MODELS],
                        "downloaded": False},
            "note": "Readability only; no model inference, transcript compatibility or shared-VEP qualification."}


def check_cadd_environments(folder, environment, names=None, include_vep=True, native=False):
    records = {}
    selected = CADD_ENVIRONMENTS if names is None else names
    for name, (prefix, major, modules) in CADD_ENVIRONMENTS.items():
        if name not in selected:
            continue
        if native:
            prefix = str(folder / "envs" / name)
            if name == "cadd-python2":
                modules = (*modules, "vcf", "statistics")
        # Blocking sockets prevents package checks from fetching missing models.
        program = """from __future__ import print_function
import hashlib, importlib, json, os, socket, sys
def no_network(*args, **kwargs):
    raise RuntimeError('Network access is disabled during CADD readiness checks')
socket.socket.connect = no_network
socket.create_connection = no_network
prefix, major, names = json.loads(sys.argv[1])
assert sys.executable == prefix + '/bin/python', 'Wrong CADD interpreter route'
assert sys.version_info[0] == major, 'Wrong CADD Python generation'
assert os.environ.get('CONDA_PREFIX') == prefix, 'Wrong CONDA_PREFIX'
assert sys.flags.no_user_site and sys.dont_write_bytecode, 'Python isolation flags missing'
versions = {}
for name in names:
    module = importlib.import_module(name)
    assert os.path.realpath(module.__file__).startswith(prefix + '/'), 'Module escaped its named environment'
    versions[name] = str(getattr(module, '__version__', 'not_exposed'))
record = {'executable': sys.executable, 'python': sys.version, 'conda_prefix': prefix,
          'versions': versions, 'network_access': 'disabled', 'model_inference': False}
if 'mmsplice' in names:
    from mmsplice.vcf_dataloader import SplicingVCFDataloader
    from mmsplice import MMSplice
    from mmsplice.utils import logit, predict_deltaLogitPsi, predict_pathogenicity, predict_splicing_efficiency
    from concise.preprocessing import encodeDNA
    root = os.path.join(os.path.dirname(importlib.import_module('mmsplice').__file__), 'models')
    assets = []
    for filename in ('Acceptor.h5', 'Donor.h5', 'Exon.h5', 'Exon_prime3.h5', 'Intron3.h5', 'Intron5.h5',
                     'Pathogenicity.pkl', 'linear_model.pkl', 'splicing_efficiency.pkl', 'grch38_exons.csv.gz'):
        path = os.path.join(root, filename)
        assert os.path.isfile(path) and os.path.getsize(path) > 0, 'Missing MMSplice packaged asset: ' + path
        with open(path, 'rb') as stream:
            stream.read(64)
        assets.append({'path': path, 'bytes': os.path.getsize(path), 'check': 'readable'})
    record['packaged_models'] = assets
print(json.dumps(record, sort_keys=True))
"""
        command = [str(folder / "bin" / name), "-c", program, json.dumps([prefix, major, modules])]
        print(f"Checking named CADD environment: {name}", flush=True)
        try:
            completed = subprocess.run(command, check=True, capture_output=True, text=True,
                                       timeout=180, env=environment)
        except subprocess.CalledProcessError as error:
            print((error.stdout or "") + (error.stderr or ""), flush=True)
            raise RuntimeError(f"Named CADD environment failed its offline import/route checks: {name}") from error
        if completed.stderr:
            print(completed.stderr, end="", flush=True)
        record = json.loads(completed.stdout.strip().splitlines()[-1])
        if record.get("executable") != prefix + "/bin/python" or record.get("conda_prefix") != prefix:
            raise RuntimeError(f"Named CADD environment returned a different route: {name}")
        records[name] = record
    if not include_vep:
        return records
    # VEP110 rejects --version; its supported --help prints the version banner
    # and exits without reading or annotating any variants.
    try:
        completed = subprocess.run([str(folder / "bin/cadd-vep"), "--help"], check=True,
                                   capture_output=True, text=True, timeout=60, env=environment)
    except subprocess.CalledProcessError as error:
        print((error.stdout or "") + (error.stderr or ""), flush=True)
        raise RuntimeError("Packaged CADD VEP failed its reference-only help/version check") from error
    version = (completed.stdout + completed.stderr).strip()
    match = re.search(r"^\s*ensembl-vep\s*:\s*(110(?:\.\d+)?)\s*$", version, re.MULTILINE)
    if match is None:
        raise RuntimeError("Packaged CADD VEP is not the expected reference-only VEP110")
    records["cadd-vep"] = {"version": match.group(1), "version_output": version,
                           "role": "reference_only_not_shared_VEP113"}
    return records

