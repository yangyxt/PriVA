#!/usr/bin/env python3
"""Install shared VEP plugin inputs while reporting gated resources explicitly."""
import argparse
import fcntl
import gzip
import json
import os
from pathlib import Path
import shutil
import subprocess
import tempfile
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
    with (args.state_dir / "shared-plugin-inputs.install.lock").open("w") as lock:
        fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        receipt = {"status": "inputs_deployed_not_pipeline_qualified", "manifest_sha256": checksum(args.manifest),
                   "files": [], "repositories": manifest["repositories"], "pending": manifest["pending"],
                   "slurm_job_id": os.environ.get("SLURM_JOB_ID")}
        for repo in manifest["repositories"]:
            target = args.cache_root / repo["path"]
            if not target.exists():
                target.parent.mkdir(parents=True, exist_ok=True)
                run(["git", "clone", "--no-checkout", repo["url"], target])
                run(["git", "-C", target, "checkout", "--detach", repo["commit"]])
            actual = subprocess.check_output(["git", "-C", str(target), "rev-parse", "HEAD"], text=True).strip()
            if actual != repo["commit"]:
                raise RuntimeError(f"Existing plugin repository has a different revision: {target}")
            run(["git", "-C", target, "diff", "--quiet", "HEAD"])
        for entry in manifest["files"]:
            target = args.cache_root / entry["path"]
            identity = target.with_name(target.name + ".origin.json")
            if target.exists() and not entry.get("md5"):
                if not identity.exists():
                    raise RuntimeError(f"Existing resource lacks a verified local identity: {target}")
                previous = json.loads(identity.read_text())
                if previous.get("url") != entry["url"] or previous.get("sha256") != checksum(target):
                    raise RuntimeError(f"Existing resource changed: {target}")
            download(entry["url"], target, entry.get("md5"), entry["bytes"])
            item = {"path": str(target), "url": entry["url"], "sha256": checksum(target), "publisher_md5": entry.get("md5")}
            write_json(identity, item)
            if "tabix" in entry and not Path(str(target) + ".tbi").exists():
                run([args.tabix, *entry["tabix"], target])
            if "gunzip_to" in entry:
                expanded = args.cache_root / entry["gunzip_to"]
                temporary = expanded.with_name(expanded.name + ".partial")
                with gzip.open(target, "rb") as src, temporary.open("wb") as dst:
                    shutil.copyfileobj(src, dst, 8 * 1024 * 1024)
                if expanded.exists() and checksum(expanded) != checksum(temporary):
                    raise RuntimeError(f"Existing expanded resource differs: {expanded}")
                temporary.replace(expanded)
            receipt["files"].append(item)
        write_json(args.state_dir / "shared-plugin-inputs.ready.json", receipt)
        print(json.dumps(receipt), flush=True)


if __name__ == "__main__":
    main()
