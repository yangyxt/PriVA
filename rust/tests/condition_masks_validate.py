#!/usr/bin/env python3
"""Freeze and validate the optional condition-mask boundary inside Slurm.

prepare --repo REPO --work-root DATA/validation/TAG --build-root DATA/build/TAG
execute --work-root DATA/validation/TAG --pins-sha256 SHA
This never installs/selects a backend, edits reference science, or touches installer caches.
"""
import argparse
import datetime
import hashlib
import json
import os
import shutil
import subprocess
import sys
from pathlib import Path

Y = Path("/stornext/snfs170/next-gen/scratch/yangyxt")
D = Path("/stornext/snfs170/next-gen/scratch/yangyxt_dt/projects/PriVA_upgrade")
TOOLCHAIN = Path("/users/u266734/.rustup/toolchains/1.90.0-x86_64-unknown-linux-gnu/bin")
PYTHON = Y / "mamba/envs/priva_acmg/bin/python"
REGISTRY = Y / "softwares/rust/cargo/registry"


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def identity(path):
    stat = path.stat()
    return stat.st_dev, stat.st_ino, stat.st_size, stat.st_mtime_ns, stat.st_ctime_ns


def utc():
    return datetime.datetime.now(datetime.timezone.utc).isoformat()


def slurm():
    assert os.environ.get("SLURM_JOB_ID") and os.environ.get("SLURM_STEP_ID"), "run inside a real Slurm step"
    return {"job_id": os.environ["SLURM_JOB_ID"], "step_id": os.environ["SLURM_STEP_ID"]}


def verify(root, expected):
    path = root / "source-pins.json"
    assert digest(path) == expected, "source manifest pin changed"
    pins = json.loads(path.read_text())
    for name, sha in pins["source_files"].items():
        assert digest(root / "source" / name) == sha, name
    for name, sha in pins["toolchain"].items():
        assert digest(Path(name)) == sha, name
    return pins


def invoke(command, log, env=None):
    with log.open("x") as stream:
        result = subprocess.run(list(map(str, command)), stdout=stream, stderr=subprocess.STDOUT, env=env)
    assert result.returncode == 0, f"command failed ({result.returncode}); retained {log}"


def prepare(args):
    native = slurm()
    assert args.repo.resolve().is_relative_to(Y)
    for path in [args.work_root, args.build_root]:
        assert path.resolve().is_relative_to(D) and not path.exists(), "fresh contained data paths required"
    args.work_root.mkdir(parents=True)
    args.build_root.mkdir(parents=True)
    source = args.work_root / "source"
    files = [args.repo / "rust" / name for name in ["Cargo.toml", "Cargo.lock", "rust-toolchain.toml", "source_functions.tsv"]]
    files += sorted((args.repo / "rust/src").rglob("*.rs"))
    files += sorted((args.repo / "rust/tests").glob("*.rs"))
    files += sorted((args.repo / "rust/tests").glob("*.py"))
    files += [args.repo / "rust/src/acmg_variant_mechanism/PANDAS_LICENSE"]
    files += [args.repo / "scripts" / name for name in ["acmg_variant_mechanism.py", "acmg_scoring.py"]]
    pins = {"created_utc": utc(), "prepare_native": native, "repo_path": str(args.repo),
            "build_root": str(args.build_root), "source_files": {}, "toolchain": {},
            "scope": "independent optional Rust boundary; main-install job 1524365 uses its separate frozen 5:09 PM CT source and does not build/select this path",
            "production_backend_selected": False, "scientific_acceptance": False}
    for path in files:
        assert not path.is_symlink(), f"source alias not allowed: {path}"
        before = identity(path)
        sha = digest(path)
        name = str(path.relative_to(args.repo))
        target = source / name
        target.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(path, target)
        assert digest(target) == sha and identity(path) == before and digest(path) == sha, name
        pins["source_files"][name] = sha
    for name in ["cargo", "rustc", "rustdoc"]:
        path = TOOLCHAIN / name
        pins["toolchain"][str(path)] = digest(path)
    assert pins["source_files"]["scripts/acmg_variant_mechanism.py"] == "c7e8da80ca0ec47fb38a5fd1d965641c1791189cb0d25a651253f44508d0bf26"
    (args.work_root / "source-pins.json").write_text(json.dumps(pins, indent=2) + "\n")
    expected = digest(args.work_root / "source-pins.json")
    verify(args.work_root, expected)
    logs = args.work_root / "logs"
    logs.mkdir()
    invoke([PYTHON, "-B", source / "rust/tests/condition_masks_oracle.py", "capture", "--source-file",
            source / "scripts/acmg_variant_mechanism.py", "--output", args.work_root / "oracle"], logs / "condition-oracle.log")
    print(json.dumps({"pins_sha256": expected, "source_files": len(pins["source_files"]),
        "oracle_manifest_sha256": digest(args.work_root / "oracle/manifest.json"), "native": native}))


def execute(args):
    native = slurm()
    pins = verify(args.work_root, args.pins_sha256)
    build = Path(pins["build_root"])
    logs = args.work_root / "logs"
    source = args.work_root / "source"
    cargo_home = build / "cargo-home"
    cargo_home.mkdir()
    (cargo_home / "registry").symlink_to(REGISTRY)
    temp = build / "tmp"
    temp.mkdir()
    env = dict(os.environ, CARGO_HOME=str(cargo_home), CARGO_TARGET_DIR=str(build / "target"),
               CARGO_BUILD_JOBS="2", RUSTC=str(TOOLCHAIN / "rustc"), RUSTDOC=str(TOOLCHAIN / "rustdoc"),
               PATH=str(TOOLCHAIN) + ":" + str(PYTHON.parent) + ":/usr/bin:/bin",
               POLARS_MAX_THREADS="1", OMP_NUM_THREADS="1", OPENBLAS_NUM_THREADS="1",
               PYTHONDONTWRITEBYTECODE="1", TMPDIR=str(temp), RUSTFLAGS="-Dwarnings")
    report = {"started_utc": utc(), "native": native, "source_pins_sha256": args.pins_sha256,
              "completed": False, "production_backend_selected": False, "scientific_acceptance": False,
              "scope": pins["scope"]}
    try:
        invoke([TOOLCHAIN / "rustc", "--version", "--verbose"], logs / "toolchain.log", env)
        invoke([TOOLCHAIN / "cargo", "test", "--locked", "--offline", "--manifest-path",
                source / "rust/Cargo.toml", "-j", "2"], logs / "cargo-tests.log", env)
        binary = args.work_root / "priva"
        shutil.copy2(build / "target/debug/priva", binary)
        report["binary_sha256"] = digest(binary)
        oracle_manifest_sha = digest(args.work_root / "oracle/manifest.json")
        invoke([PYTHON, "-B", source / "rust/tests/condition_masks_compare.py", "compare", "--binary", binary,
                "--binary-sha256", report["binary_sha256"], "--oracle", args.work_root / "oracle",
                "--oracle-manifest-sha256", oracle_manifest_sha, "--output", args.work_root / "comparison"],
               logs / "condition-comparison.log", env)
        invoke([PYTHON, "-B", source / "rust/tests/mechanism_masks_oracle.py", "capture", "--output",
                args.work_root / "mechanism-oracle"], logs / "mechanism-oracle.log", env)
        mechanism = args.work_root / "mechanism-oracle"
        manifest = {p.name: digest(p) for p in sorted(mechanism.iterdir()) if p.is_file()}
        manifest_path = args.work_root / "mechanism-oracle-manifest.json"
        manifest_path.write_text(json.dumps(manifest, indent=2) + "\n")
        invoke([PYTHON, "-B", source / "rust/tests/mechanism_masks_compare.py", "compare", "--binary", binary,
                "--oracle", mechanism, "--oracle-manifest", manifest_path, "--output", args.work_root / "mechanism-comparison"],
               logs / "mechanism-comparison.log", env)
        invoke([PYTHON, "-B", source / "rust/tests/differential_scoring.py", binary], logs / "scoring-comparison.log", env)
        verify(args.work_root, args.pins_sha256)
        assert digest(binary) == report["binary_sha256"]
        report["completed"] = True
    except Exception as error:
        report["failure"] = {"type": type(error).__name__, "message": str(error)}
        raise
    finally:
        report["finished_utc"] = utc()
        with (args.work_root / "validation.json").open("x") as stream:
            json.dump(report, stream, indent=2)
            stream.write("\n")
        print(json.dumps(report), flush=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest="command", required=True)
    prepare_parser = sub.add_parser("prepare")
    for name in ["repo", "work-root", "build-root"]:
        prepare_parser.add_argument("--" + name, type=Path, required=True)
    execute_parser = sub.add_parser("execute")
    execute_parser.add_argument("--work-root", type=Path, required=True)
    execute_parser.add_argument("--pins-sha256", required=True)
    args = parser.parse_args()
    (prepare if args.command == "prepare" else execute)(args)


if __name__ == "__main__":
    main()
