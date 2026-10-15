#!/usr/bin/env python3
"""Resumable, separately locked cluster resource stages for the GRCh37 baseline."""
import argparse
import fcntl
import gzip
import hashlib
import json
import os
from pathlib import Path
import shutil
import subprocess
import tempfile
from urllib.parse import urlparse

REPO = Path(__file__).resolve().parents[1]


def run(command, **kwargs):
    print("+ " + " ".join(map(str, command)), flush=True)
    subprocess.run(list(map(str, command)), check=True, **kwargs)


def checksum(path, algorithm="sha256"):
    digest = hashlib.new(algorithm)
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def download(url, destination, md5=None, expected_size=None):
    """Never publish partial downloads; recheck existing artifacts before reuse."""
    destination.parent.mkdir(parents=True, exist_ok=True)
    partial = destination.with_name(destination.name + ".partial")
    candidate = destination if destination.exists() else partial
    if candidate == partial:
        run(["curl", "--fail", "--location", "--retry", "3", "--connect-timeout", "60",
             "--continue-at", "-", "--output", partial, url])
    if md5 and checksum(candidate, "md5") != md5:
        raise RuntimeError(f"Checksum mismatch: {candidate}; retained for investigation")
    if expected_size and candidate.stat().st_size != expected_size:
        raise RuntimeError(f"Unexpected size: {candidate}")
    if candidate == partial:
        partial.replace(destination)
    return destination


def write_json(path, value):
    temporary = path.with_name(path.name + ".partial")
    temporary.write_text(json.dumps(value, indent=2) + "\n")
    temporary.replace(path)


def source(args, manifest):
    destination = args.cache_root / "CADD/CADD-scripts-1.7.2"
    destination.parent.mkdir(parents=True, exist_ok=True)
    if not destination.exists():
        run(["git", "clone", "--no-checkout", manifest["cadd"]["source_url"], destination])
        run(["git", "-C", destination, "checkout", "--detach", manifest["cadd"]["source_commit"]])
    actual = subprocess.check_output(["git", "-C", str(destination), "rev-parse", "HEAD"], text=True).strip()
    if actual != manifest["cadd"]["source_commit"]:
        raise RuntimeError("CADD source revision does not match pinned manifest")
    run(["git", "-C", destination, "diff", "--quiet", "HEAD", "--", "src", "config", "envs", "Snakefile"])
    return destination


def reference(args, manifest):
    folder = args.cache_root / "references"
    compressed = download(manifest["reference"]["url"], folder / "ucsc.hg19.fasta.gz", manifest["reference"]["md5"])
    target = folder / "ucsc.hg19.fasta"
    if not target.exists():
        partial = target.with_name(target.name + ".partial")
        with gzip.open(compressed, "rb") as src, partial.open("wb") as dest:
            shutil.copyfileobj(src, dest, 8 * 1024 * 1024)
        partial.replace(target)
    if not Path(str(target) + ".fai").exists():
        run([args.samtools, "faidx", target])
    return {"reference": str(target), "sha256": checksum(target), "archive_md5": manifest["reference"]["md5"]}


def vep(args, manifest):
    folder = args.cache_root / "VEP_cache"
    version = manifest["vep"]["cache_release"]
    target = folder / f"homo_sapiens_merged/{version}_GRCh37"
    archive = download(manifest["vep"]["url"], args.cache_root / ".downloads" / Path(urlparse(manifest["vep"]["url"]).path).name,
                       expected_size=manifest["vep"]["bytes"])
    published = manifest["vep"]["bsd_sum"]
    measured = subprocess.check_output(["sum", "-r", str(archive)], text=True).split()[:2]
    if list(map(int, measured)) != published:
        raise RuntimeError("VEP archive differs from publisher CHECKSUMS")
    if target.exists() and not (target / "info.txt").exists():
        raise RuntimeError(f"Incomplete existing VEP cache: {target}")
    if not target.exists():
        folder.mkdir(parents=True, exist_ok=True)
        staging = Path(tempfile.mkdtemp(prefix=".extract-", dir=folder))
        run(["tar", "-xzf", archive, "-C", staging])
        extracted = staging / "homo_sapiens_merged" / f"{version}_GRCh37"
        if not (extracted / "info.txt").exists():
            raise RuntimeError(f"VEP archive missing expected cache: {extracted}")
        target.parent.mkdir(parents=True, exist_ok=True)
        extracted.rename(target)
        shutil.rmtree(staging)
    return {"cache": str(target), "archive_sha256": checksum(archive), "assembly": "GRCh37", "release": version}


def cadd_runtime(args, manifest):
    location = source(args, manifest)
    folder = args.software_root / "cadd-1.7.2"
    folder.mkdir(parents=True, exist_ok=True)
    image = folder / "cadd-upstream-0.1.1.sif"
    environment = dict(os.environ)
    environment["APPTAINER_CACHEDIR"] = str(folder / "image-cache")
    environment["APPTAINER_TMPDIR"] = str(args.tmp_root / "cadd-apptainer")
    Path(environment["APPTAINER_TMPDIR"]).mkdir(parents=True, exist_ok=True)
    if not image.exists():
        partial = image.with_name(image.name + ".partial")
        run(["apptainer", "pull", "--force", partial, manifest["cadd"]["image"]], env=environment)
        partial.replace(image)
    run(["apptainer", "exec", "--cleanenv", image,
         "/conda-envs/fe7049a03d989a2faa6d7688f5987d09/bin/python", "-c",
         "import sys, sklearn, numpy, pandas; print(sys.version); print(sklearn.__version__)"], env=environment)
    return {"source": str(location), "image": str(image), "image_sha256": checksum(image), "image_source": manifest["cadd"]["image"]}


def cadd_annotations(args, manifest):
    location = source(args, manifest)
    folder = location / "data/annotations"
    folder.mkdir(parents=True, exist_ok=True)
    target = folder / "GRCh37_v1.7"
    archive = download(manifest["cadd"]["annotation_url"], args.cache_root / ".downloads/GRCh37_v1.7.tar.gz",
                       manifest["cadd"]["annotation_md5"])
    if target.exists() and not (target / "vep").is_dir():
        raise RuntimeError(f"Incomplete existing CADD annotations: {target}")
    if not target.exists():
        staging = Path(tempfile.mkdtemp(prefix=".extract-", dir=folder))
        run(["tar", "-xzf", archive, "-C", staging])
        extracted = staging / "GRCh37_v1.7"
        if not (extracted / "vep").is_dir():
            raise RuntimeError("CADD annotation archive lacks the expected VEP resources")
        extracted.rename(target)
        shutil.rmtree(staging)
    return {"annotations": str(target), "archive_md5": manifest["cadd"]["annotation_md5"]}


def cadd_prescores(args, manifest):
    location = source(args, manifest)
    folder = location / "data/prescored/GRCh37_v1.7/incl_anno"
    for filename, md5 in manifest["cadd"]["prescores"].items():
        download(manifest["cadd"]["prescore_url_root"] + filename, folder / filename, md5)
    return {"prescores": str(folder), "md5": manifest["cadd"]["prescores"], "note": manifest["cadd"]["checksum_note"]}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("stage", choices=["reference", "vep", "cadd-runtime", "cadd-annotations", "cadd-prescores"])
    parser.add_argument("--manifest", type=Path, default=REPO / "envs/cluster-resources.hg19.json")
    parser.add_argument("--cache-root", type=Path, required=True)
    parser.add_argument("--software-root", type=Path, required=True)
    parser.add_argument("--state-dir", type=Path, required=True)
    parser.add_argument("--tmp-root", type=Path, required=True)
    parser.add_argument("--samtools", default="samtools")
    args = parser.parse_args()
    manifest = json.loads(args.manifest.read_text())
    if manifest["assembly"] != "GRCh37":
        parser.error("This reviewed installer currently supports the GRCh37 baseline only")
    args.state_dir.mkdir(parents=True, exist_ok=True)
    with (args.state_dir / (args.stage + ".install.lock")).open("w") as lock:
        fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        # Source creation is shared by CADD stages; serialize only this step.
        if args.stage.startswith("cadd-"):
            with (args.state_dir / "cadd-source.install.lock").open("w") as source_lock:
                fcntl.flock(source_lock, fcntl.LOCK_EX)
                source(args, manifest)
        function = globals()[args.stage.replace("-", "_")]
        result = function(args, manifest)
        result.update({"stage": args.stage, "status": "deployed", "slurm_job_id": os.environ.get("SLURM_JOB_ID"),
                       "manifest_sha256": checksum(args.manifest)})
        write_json(args.state_dir / (args.stage + ".ready.json"), result)
        print(json.dumps(result), flush=True)


if __name__ == "__main__":
    main()
