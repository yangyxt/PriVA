#!/usr/bin/env python3
"""Download the pinned gnomAD source set; no implicit hg19 conversion."""
import argparse
import fcntl
import json
import os
from pathlib import Path
from deploy_cluster_resources import checksum, download, run, write_json


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest", type=Path, required=True)
    parser.add_argument("--cache-root", type=Path, required=True)
    parser.add_argument("--state-dir", type=Path, required=True)
    parser.add_argument("--tabix", type=Path, required=True)
    args = parser.parse_args()
    manifest = json.loads(args.manifest.read_text())
    args.state_dir.mkdir(parents=True, exist_ok=True)
    with (args.state_dir / "gnomad-source.install.lock").open("w") as lock:
        fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        for entry in manifest["files"]:
            target = args.cache_root / "gnomAD" / entry["path"]
            download(entry["url"], target, entry["md5"], entry["bytes"])
            if entry["path"].endswith(".tbi"):
                run([args.tabix, "-l", str(target)[:-4]])
        write_json(args.state_dir / "gnomad-source.ready.json", {
            "status": "source_inputs_deployed", "assembly": "GRCh38", "release": "4.1",
            "scope": "GRCh37 liftover and calibration caches are not yet qualified",
            "manifest_sha256": checksum(args.manifest), "files": manifest["files"],
            "slurm_job_id": os.environ.get("SLURM_JOB_ID")})


if __name__ == "__main__":
    main()
