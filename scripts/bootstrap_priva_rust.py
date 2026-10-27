#!/usr/bin/env python3
"""Install/check the optional PriVA procedure CLI without selecting a backend.

  python scripts/bootstrap_priva_rust.py install --source-root REPO --cache-root CACHE
  python scripts/bootstrap_priva_rust.py check --source-root REPO --cache-root CACHE

An explicit qualified import additionally requires --qualified-binary,
--build-pins/--build-pins-sha256 and --validation/--validation-sha256. Import
retains the original executable's actual profile and scoped qualification;
it never turns software tests into family or clinical acceptance. A cold
build uses Cargo.lock and the exact rust-toolchain.toml release, with separate
provenance and no inherited scientific qualification. All build/download/log
products stay below cache_root/runtime/priva-rust. No environment, shell,
PATH, ranking configuration or original PriVA file is changed.
"""
from __future__ import annotations

import argparse
import contextlib
from datetime import datetime, timezone
import fcntl
import hashlib
import json
import os
from pathlib import Path
import platform
import re
import shutil
import stat
import subprocess
import sys
import tempfile
import tomllib
import urllib.request

KIND = "priva-managed-rust-cli-v1"
HOST = "x86_64-unknown-linux-gnu"
RUSTUP_URL = "https://static.rust-lang.org/rustup/archive/1.28.2/" + HOST + "/rustup-init"
RUSTUP_SHA256 = "20a06e644b0d9bd2fbdbfd52d42540bdde820ea7df86e92e533c073da0cdd43c"
IMPORTED_PROFILE = "test/dev, opt-level=0, debug=0; no fast-math"
SOURCE_FIXED = ("rust/Cargo.toml", "rust/Cargo.lock", "rust/rust-toolchain.toml")


def require(ok, message):
    if not ok:
        raise ValueError(message)


def canonical(value):
    return json.dumps(value, sort_keys=True, separators=(",", ":"), allow_nan=False).encode()


def sha(raw):
    return hashlib.sha256(raw).hexdigest()


def identity(path):
    """Cross-node persisted observations deliberately omit node-local st_dev."""
    s = Path(path).lstat()
    return dict(size=s.st_size, inode=s.st_ino, mode=s.st_mode, nlink=s.st_nlink,
                mtime_ns=s.st_mtime_ns, ctime_ns=s.st_ctime_ns)


def fd_identity(fd):
    s = os.fstat(fd)
    return dict(size=s.st_size, inode=s.st_ino, mode=s.st_mode, nlink=s.st_nlink,
                mtime_ns=s.st_mtime_ns, ctime_ns=s.st_ctime_ns)


def no_links(path, root=None):
    path = Path(path)
    require(path.is_absolute() and ".." not in path.parts, "Absolute non-traversing path required")
    root = Path(path.anchor) if root is None else Path(root)
    require(path.is_relative_to(root), "Path escapes its managed root")
    current = path
    while current != root:
        require(not current.is_symlink(), "Refusing symlinked path component: " + str(current))
        current = current.parent


def read_stable(path, limit=None):
    path = Path(path); no_links(path)
    before = identity(path)
    require(stat.S_ISREG(before["mode"]), "Expected ordinary file: " + str(path))
    require(limit is None or before["size"] <= limit, "Unexpectedly large metadata")
    with os.fdopen(os.open(path, os.O_RDONLY | os.O_NOFOLLOW), "rb") as stream:
        require(fd_identity(stream.fileno()) == before, "Input changed before open")
        raw = stream.read()
    require(identity(path) == before and len(raw) == before["size"], "Input changed during read")
    return raw, before


def digest(path):
    path = Path(path); no_links(path)
    before = identity(path); require(stat.S_ISREG(before["mode"]), "Expected ordinary content")
    h = hashlib.sha256()
    with os.fdopen(os.open(path, os.O_RDONLY | os.O_NOFOLLOW), "rb") as stream:
        require(fd_identity(stream.fileno()) == before, "Content changed before open")
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            h.update(block)
    require(identity(path) == before, "Content changed during hashing")
    return h.hexdigest()


def load(path):
    raw, before = read_stable(path, 16 * 1024**2)
    return json.loads(raw), raw, before


def save_new(path, value):
    with Path(path).open("xb") as stream:
        stream.write(canonical(value) + b"\n"); stream.flush(); os.fsync(stream.fileno())


def atomic_new(path, raw):
    """Publish a complete ordinary file without replacement; SNFS supports link.

    renameat2(RENAME_NOREPLACE) returns EINVAL on the deployed SNFS. Hardlinks
    here concern only NEW owned metadata/staging copies, never an import input.
    """
    path = Path(path); no_links(path.parent)
    fd, temporary = tempfile.mkstemp(prefix=".publish-", dir=path.parent)
    try:
        with os.fdopen(fd, "wb") as stream:
            stream.write(raw); stream.flush(); os.fsync(stream.fileno())
        os.link(temporary, path, follow_symlinks=False)
    finally:
        os.unlink(temporary)


def publication_path(target):
    return target.parent / ("publication-" + target.name + ".json")


def target_creation(target):
    no_links(target)
    observed = target.lstat()
    require(stat.S_ISDIR(observed.st_mode), "Publication target is not an ordinary directory")
    # Contents change directory times/link counts; inode, owner and mode bind
    # the directory created by this transaction across its member writes.
    return dict(inode=observed.st_ino, uid=observed.st_uid, mode=observed.st_mode)


def publication_ledger(target, require_created=False):
    ledger, raw, _ = load(publication_path(target))
    legacy = {"kind", "identity", "stage", "stage_inode", "members"}
    require(set(ledger) in (legacy, legacy | {"protocol_version"})
            and ledger["kind"] == KIND + "-publication" and ledger["identity"] == target.name,
            "Foreign publication intent")
    if "protocol_version" not in ledger:
        # Completed earlier artifacts remain byte-bound and read-only reusable;
        # an old incomplete intent cannot prove who created its target.
        require(not require_created, "Incomplete legacy publication has no target-creation proof")
    else:
        require(ledger["protocol_version"] == 2, "Unsupported publication protocol")
        created_path = target.parent / ("publication-" + target.name + ".created.json")
        require(created_path.exists(), "Publication target creation was not confirmed; retained")
        created, _, _ = load(created_path)
        require(created == dict(kind=KIND + "-target-created", identity=target.name,
                                publication_sha256=sha(raw), target=target_creation(target)),
                "Publication target-creation proof differs; retained")
    return ledger


def resume_publication(target):
    """Resume only an exact owned member ledger; commit readiness last.

    A directory or partial file is never sufficient for reuse. The ready
    receipt is linked atomically only after every other declared byte/mode
    matches. Existing unowned paths, extra members and links are refused.
    """
    target = Path(target); base = target.parent; no_links(target)
    ledger = publication_ledger(target, require_created=True)
    stage = base / ledger["stage"]
    require(stage.parent == base and stage.name.startswith(".artifact-work-") and stage.is_dir()
            and not stage.is_symlink() and stage.stat().st_ino == ledger["stage_inode"],
            "Owned prepared publication is unavailable")
    members = ledger["members"]
    require(set(members) in ({"bin/priva", "receipt.json", "evidence/build-pins.json", "evidence/validation.json"},
                            {"bin/priva", "receipt.json", "evidence/build.json", "evidence/build.log"}),
            "Unexpected publication member set")
    expected = set(members) | {"bin", "evidence"}
    require(target.is_dir() and not target.is_symlink()
            and {str(path.relative_to(target)) for path in target.rglob("*")} <= expected,
            "Foreign or extra partial publication members")
    for folder in ("bin", "evidence"):
        path = target / folder; no_links(path, base); path.mkdir(exist_ok=True)
        require(path.is_dir(), "Partial publication directory differs")
    # A ready receipt is the transaction's commit, never a partially written file.
    order = sorted(set(members) - {"receipt.json"}) + ["receipt.json"]
    for name in order:
        source, destination = stage / name, target / name
        row = members[name]; no_links(source, base); no_links(destination, base)
        require(identity(source)["mode"] == row["mode"] and digest(source) == row["sha256"],
                "Prepared publication bytes/mode changed")
        try:
            os.link(source, destination, follow_symlinks=False)
        except FileExistsError:
            require(identity(destination)["mode"] == row["mode"] and digest(destination) == row["sha256"],
                    "Existing partial publication differs; retained")
    # Successful fresh staging copies can be removed; a failed attempt's stage
    # and failure record remain intact for review and recovery provenance.
    if not (stage / "failure.json").exists():
        require({str(path.relative_to(stage)) for path in stage.rglob("*")}
                == expected | {"publication.json"}, "Unexpected staging member; retained")
        for name in (*members, "publication.json"):
            (stage / name).unlink()
        (stage / "bin").rmdir(); (stage / "evidence").rmdir(); stage.rmdir()


def publish_artifact(stage, target):
    members = {str(path.relative_to(stage)): dict(sha256=digest(path), mode=identity(path)["mode"])
               for path in stage.rglob("*") if path.is_file() and not path.is_symlink()}
    ledger = dict(kind=KIND + "-publication", protocol_version=2, identity=target.name,
                  stage=stage.name, stage_inode=stage.stat().st_ino, members=members)
    raw = canonical(ledger) + b"\n"
    (stage / "publication.json").write_bytes(raw)
    # An intent alone never owns an empty racing directory. Persist positive
    # proof only AFTER our exclusive mkdir succeeds and BEFORE member writes.
    # An interruption in this small gap is retained and fails closed on retry.
    atomic_new(publication_path(target), raw)
    target.mkdir(mode=0o700)  # Refuses an independently created racing target.
    atomic_new(target.parent / ("publication-" + target.name + ".created.json"),
               canonical(dict(kind=KIND + "-target-created", identity=target.name,
                              publication_sha256=sha(raw), target=target_creation(target))) + b"\n")
    resume_publication(target)


def source_state(root):
    root = Path(root); no_links(root)
    require((root / "rust/src/main.rs").is_file(), "Expected PriVA Rust source root")
    paths = [root / name for name in SOURCE_FIXED]
    # Keep non-Rust source companions such as PANDAS_LICENSE in the same
    # producer identity; a filename suffix is not a complete source manifest.
    paths += sorted(path for path in (root / "rust/src").rglob("*")
                    if path.is_file() or path.is_symlink())
    if (root / "rust/build.rs").exists():
        paths.append(root / "rust/build.rs")
    # A custom Cargo configuration is a build input, not an invisible override.
    for config in (root / "rust/.cargo/config", root / "rust/.cargo/config.toml"):
        require(not config.exists() and not config.is_symlink(), "Custom Cargo configuration needs explicit build support")
    files = {str(path.relative_to(root)): digest(path) for path in paths}
    manifest = tomllib.loads(read_stable(root / "rust/Cargo.toml")[0].decode())
    toolchain = tomllib.loads(read_stable(root / "rust/rust-toolchain.toml")[0].decode())["toolchain"]
    channel = toolchain["channel"]
    require(re.fullmatch(r"\d+\.\d+\.\d+", channel) is not None, "Pin an exact Rust release")
    require(manifest["package"]["name"] == "priva-tabular"
            and any(b.get("name") == "priva" for b in manifest.get("bin", [])), "Expected PriVA CLI package")
    return dict(files=files, source_sha256=sha(canonical(files)), channel=channel,
                package_version=manifest["package"]["version"], cargo_profiles=manifest.get("profile", {}))


def managed_base(cache_root, create=False):
    require(Path(cache_root).is_absolute(), "An absolute cache_root is required")
    supplied = Path(os.path.abspath(cache_root))
    require(supplied.is_absolute() and supplied.is_dir(), "An existing cache_root is required")
    # Match the installer's established root-alias normalization, then reject
    # links in every managed child. No other root is resolved implicitly.
    cache = supplied.resolve(strict=True)
    require(cache != Path("/"), "cache_root must be a dedicated directory")
    runtime, base = cache / "runtime", cache / "runtime/priva-rust"
    no_links(base, cache)
    anchor = runtime / ".priva-rust-owner.json"
    owner_value = dict(kind=KIND, schema_version=1, cache_root=str(cache))
    if create:
        runtime.mkdir(exist_ok=True)
        no_links(runtime, cache);require(runtime.is_dir(), "Managed runtime parent is not a directory")
        require(not base.exists() or anchor.exists(), "Foreign runtime directory; retained")
        if not anchor.exists():
            try: atomic_new(anchor, canonical(owner_value) + b"\n")
            except FileExistsError: pass  # A concurrent installer must have the identical owner below.
        require(load(anchor)[0] == owner_value, "Foreign runtime ownership anchor")
        if not base.exists():
            try: base.mkdir(mode=0o700)
            except FileExistsError: pass
        no_links(base, cache)
        if not (base / "owner.json").exists():
            require(not list(base.iterdir()), "Unowned runtime members; retained")
            try: atomic_new(base / "owner.json", canonical(owner_value) + b"\n")
            except FileExistsError: pass
        require(load(base / "owner.json")[0] == owner_value, "Foreign runtime owner")
        if not (base / ".install.lock").exists():
            try:
                fd = os.open(base / ".install.lock", os.O_CREAT | os.O_EXCL | os.O_WRONLY | os.O_NOFOLLOW, 0o600)
                os.close(fd)
            except FileExistsError: pass
    require(base.is_dir(), "PriVA Rust runtime is not installed")
    require(load(anchor)[0] == owner_value, "Foreign runtime ownership anchor")
    owner, _, _ = load(base / "owner.json")
    require(owner == owner_value, "Foreign runtime ownership; retained")
    return base


@contextlib.contextmanager
def locked(base, exclusive):
    lock = base / ".install.lock"; no_links(lock, base)
    fd = os.open(lock, (os.O_RDWR if exclusive else os.O_RDONLY) | os.O_NOFOLLOW)
    try:
        s = os.fstat(fd)
        require(stat.S_ISREG(s.st_mode) and s.st_nlink == 1 and s.st_uid == os.geteuid(), "Foreign runtime lock")
        fcntl.flock(fd, fcntl.LOCK_EX if exclusive else fcntl.LOCK_SH)
        yield
    finally:
        os.close(fd)


def import_evidence(args, state):
    fields = (args.qualified_binary, args.build_pins, args.build_pins_sha256,
              args.validation, args.validation_sha256)
    require(all(fields) or not any(fields), "Qualified import requires all five explicit binary/evidence pins")
    if not any(fields):
        return None
    pins, raw_pins, pin_before = load(args.build_pins)
    result, raw_result, result_before = load(args.validation)
    require(sha(raw_pins) == args.build_pins_sha256 and sha(raw_result) == args.validation_sha256,
            "Explicit import receipt checksum differs")
    expected = {name: value for name, value in pins["source_files"].items()
                if name.startswith("rust/src/") or name in SOURCE_FIXED or name == "rust/build.rs"}
    require(expected == state["files"], "Qualified CLI source/Cargo.lock/toolchain source differs")
    require(result["completed"] is True and result["source_pins_sha256"] == sha(raw_pins)
            and type(result["rust_tests_passed"]) is int and result["rust_tests_passed"] > 0,
            "Completed source-bound build qualification required")
    require(pins["actual_Cargo_profile"] == IMPORTED_PROFILE
            and state["cargo_profiles"].get("dev", {}) == {"debug": 0}
            and state["cargo_profiles"].get("test", {}) == {"debug": 0},
            "Import supports only the explicitly qualified actual dev/test profile")
    tool_files = {}
    for name, wanted in pins["toolchain"].items():
        key = Path(name).name
        require(key in ("cargo", "rustc", "rustdoc", "rustfmt") and key not in tool_files
                and re.fullmatch("[0-9a-f]{64}", wanted), "Invalid qualified toolchain pins")
        tool_files[key] = wanted
    require(set(tool_files) == {"cargo", "rustc", "rustdoc", "rustfmt"}, "Complete qualified toolchain pins required")
    require(re.fullmatch("[0-9a-f]{64}", result["binary_sha256"]), "Invalid qualified executable digest")
    return dict(binary=Path(args.qualified_binary), binary_sha256=result["binary_sha256"],
        evidence={"build-pins.json": raw_pins, "validation.json": raw_result},
        input_observations={str(args.build_pins): pin_before, str(args.validation): result_before},
        profile=dict(name="dev", actual=IMPORTED_PROFILE, opt_level=0, debug=0, fast_math=False),
        toolchain=dict(channel=state["channel"], host=HOST, files=tool_files),
        qualification=dict(origin="explicit_qualified_import", native=result["native"],
            scope=result["scope"], rust_tests_passed=result["rust_tests_passed"],
            scientific_acceptance=result["scientific_acceptance"],
            production_backend_selected=result["production_backend_selected"]))


def probe(binary, version):
    output = subprocess.check_output([str(binary), "--version"], text=True, timeout=30,
                                     env={k: v for k, v in os.environ.items() if k not in ("LD_PRELOAD", "LD_LIBRARY_PATH")})
    require(output.startswith("priva " + version + " "), "Unexpected PriVA CLI version")
    return output.strip()


def artifact_descriptor(state, spec):
    return dict(source_files=state["files"], source_sha256=state["source_sha256"],
                toolchain=spec["toolchain"], profile=spec["profile"], binary_sha256=spec["binary_sha256"])


def verify_artifact(root, state):
    no_links(root)
    publication = publication_ledger(root)
    receipt, raw_receipt, before = load(root / "receipt.json")
    require(receipt["kind"] == KIND and receipt["schema_version"] == 1, "Foreign artifact receipt")
    descriptor = receipt["artifact"]
    require(sha(canonical(descriptor)) == root.name and descriptor["source_files"] == state["files"]
            and descriptor["source_sha256"] == state["source_sha256"]
            and descriptor["toolchain"]["channel"] == state["channel"]
            and descriptor["toolchain"]["host"] == HOST, "Installed source/toolchain/profile identity differs")
    evidence = receipt["evidence_sha256"]
    require(set(evidence) in ({"build-pins.json", "validation.json"}, {"build.json", "build.log"}), "Unexpected build evidence")
    for name, wanted in evidence.items():
        require(digest(root / "evidence" / name) == wanted, "Installed build evidence changed")
    if "build-pins.json" in evidence:
        pins, _, _ = load(root / "evidence/build-pins.json")
        result, _, _ = load(root / "evidence/validation.json")
        require(result["completed"] is True and result["source_pins_sha256"] == evidence["build-pins.json"]
                and result["binary_sha256"] == descriptor["binary_sha256"]
                and type(result["rust_tests_passed"]) is int and result["rust_tests_passed"] > 0
                and pins["actual_Cargo_profile"] == IMPORTED_PROFILE
                and descriptor["profile"] == dict(name="dev", actual=IMPORTED_PROFILE, opt_level=0, debug=0, fast_math=False)
                and state["cargo_profiles"].get("dev", {}) == {"debug": 0}
                and state["cargo_profiles"].get("test", {}) == {"debug": 0},
                "Installed import profile/qualification does not bind this executable")
        source_pins = {name: value for name, value in pins["source_files"].items()
            if name.startswith("rust/src/") or name in SOURCE_FIXED or name == "rust/build.rs"}
        require(source_pins == descriptor["source_files"], "Installed import producer source differs")
        pinned_tools = {Path(k).name: v for k, v in pins["toolchain"].items()}
        require(len(pins["toolchain"]) == len(pinned_tools) == 4
                and set(pinned_tools) == {"cargo", "rustc", "rustdoc", "rustfmt"}
                and pinned_tools == descriptor["toolchain"]["files"], "Installed compiler pins differ")
        require(receipt["qualification"] == dict(origin="explicit_qualified_import", native=result["native"],
            scope=result["scope"], rust_tests_passed=result["rust_tests_passed"],
            scientific_acceptance=result["scientific_acceptance"],
            production_backend_selected=result["production_backend_selected"]), "Original qualification was relabelled")
    else:
        built, _, _ = load(root / "evidence/build.json")
        require(built["artifact"] == descriptor and receipt["qualification"] ==
                dict(origin="cold_locked_build", scientific_acceptance=False, production_backend_selected=False),
                "Cold build cannot inherit imported qualification")
    binary = root / "bin/priva"
    require(identity(binary)["mode"] & 0o111 and digest(binary) == descriptor["binary_sha256"], "Installed executable differs")
    actual = {str(path.relative_to(root)) for path in root.rglob("*")}
    expected = {"receipt.json", "bin", "bin/priva", "evidence", *("evidence/" + name for name in evidence)}
    require(actual == expected and all(not path.is_symlink() for path in root.rglob("*")), "Unexpected artifact member or link")
    expected_hashes = {"receipt.json": sha(raw_receipt), "bin/priva": descriptor["binary_sha256"],
                       **{"evidence/" + name: wanted for name, wanted in evidence.items()}}
    require(set(publication["members"]) == set(expected_hashes)
            and all(publication["members"][name] == dict(sha256=wanted, mode=identity(root / name)["mode"])
                    for name, wanted in expected_hashes.items()), "Published member ledger differs")
    require(identity(root / "receipt.json") == before and (root / "receipt.json").read_bytes() == raw_receipt,
            "Installed receipt changed during verification")
    return dict(identity=root.name, binary=str(binary), receipt=str(root / "receipt.json"),
                profile=descriptor["profile"], qualification=receipt["qualification"], ready=True)


def find_existing(base, state, selected=None, recover=False):
    found = []
    identities = {path.name for path in base.iterdir() if re.fullmatch("[0-9a-f]{64}", path.name)}
    identities |= {path.name[12:-5] for path in base.glob("publication-*.json")
                   if re.fullmatch(r"publication-[0-9a-f]{64}[.]json", path.name)}
    for key in sorted(identities):
        path = base / key
        if not (path / "receipt.json").exists():
            ledger, _, _ = load(publication_path(path))
            stage = base / ledger["stage"]
            require(stage.parent == base and stage.name.startswith(".artifact-work-"), "Foreign prepared path")
            prepared, _, _ = load(stage / "receipt.json")
            if prepared.get("artifact", {}).get("source_sha256") != state["source_sha256"]:
                continue
            require(recover, "Owned publication is incomplete; run install to recover it")
            resume_publication(path)
        if re.fullmatch("[0-9a-f]{64}", path.name):
            require(path.is_dir() and not path.is_symlink(), "Foreign or damaged artifact directory")
            receipt, _, _ = load(path / "receipt.json")
            if receipt.get("artifact", {}).get("source_sha256") == state["source_sha256"]:
                if selected is None or selected == path.name:
                    found.append(verify_artifact(path, state))
    require(len(found) <= 1, "Several matching builds exist; select --identity explicitly")
    if selected:
        require(found, "Requested managed build does not match the current source")
    return found[0] if found else None


def toolchain_descriptor(bin_dir, state):
    bin_dir = Path(bin_dir); no_links(bin_dir)
    version = subprocess.check_output([str(bin_dir / "rustc"), "-vV"], text=True, timeout=30)
    require("release: " + state["channel"] + "\n" in version and "host: " + HOST + "\n" in version,
            "Compiler does not match the pinned Rust release/host")
    files = {name: digest(bin_dir / name) for name in ("cargo", "rustc", "rustdoc", "rustfmt")}
    return dict(channel=state["channel"], host=HOST, files=files)


def download_bootstrap(path):
    h = hashlib.sha256(); total = 0
    with urllib.request.urlopen(RUSTUP_URL, timeout=60) as source, path.open("xb") as target:
        for block in iter(lambda: source.read(1024 * 1024), b""):
            total += len(block);require(total <= 64 * 1024**2, "Unexpectedly large Rustup bootstrap")
            target.write(block); h.update(block)
    require(h.hexdigest() == RUSTUP_SHA256, "Pinned official rustup bootstrap checksum differs")
    path.chmod(0o700)


def obtain_toolchain(base, state, explicit=None, offline=False):
    if explicit:
        path = Path(explicit); toolchain_descriptor(path, state); return path
    candidates = []
    rustup = shutil.which("rustup")
    if rustup:
        result = subprocess.run([rustup, "which", "--toolchain", state["channel"], "rustc"],
                                capture_output=True, text=True, timeout=30)
        if result.returncode == 0:
            candidates.append(Path(result.stdout.strip()).parent)
    for home in (os.environ.get("RUSTUP_HOME"), str(Path.home() / ".rustup")):
        if home:
            candidates.append(Path(home) / "toolchains" / (state["channel"] + "-" + HOST) / "bin")
    for path in candidates:
        if (path / "rustc").is_file():
            try: toolchain_descriptor(path, state)
            except (OSError, ValueError, subprocess.SubprocessError): continue
            return path
    target = base / ("toolchain-" + state["channel"])
    no_links(target, base)
    intent = base / ("toolchain-" + state["channel"] + ".intent.json")
    owner_value = dict(kind=KIND + "-toolchain", channel=state["channel"], bootstrap_sha256=RUSTUP_SHA256)
    if (target / "owner.json").exists():
        owner, _, _ = load(target / "owner.json")
        require(owner == owner_value and load(intent)[0] == owner_value,
                "Foreign toolchain directory; retained")
        path = target / "rustup/toolchains" / (state["channel"] + "-" + HOST) / "bin"
        toolchain_descriptor(path, state); return path
    require(not offline, "Pinned toolchain is absent and offline mode was requested")
    require(not target.exists() or intent.exists(), "Foreign incomplete toolchain directory; retained")
    if not intent.exists(): atomic_new(intent, canonical(owner_value) + b"\n")
    require(load(intent)[0] == owner_value, "Foreign pending toolchain intent")
    if not target.exists(): target.mkdir(mode=0o700)
    no_links(target, base)
    require(all(path.name in ("rustup-init", "rustup", "cargo") or path.name.startswith("install-")
                for path in target.iterdir()), "Foreign pending toolchain member")
    init = target / "rustup-init"
    if init.exists(): require(digest(init) == RUSTUP_SHA256, "Existing bootstrap changed; retained")
    else: download_bootstrap(init)
    env = dict(os.environ, RUSTUP_HOME=str(target / "rustup"), CARGO_HOME=str(target / "cargo"))
    fd, log_name = tempfile.mkstemp(prefix="install-", suffix=".log", dir=target)
    with os.fdopen(fd, "wb") as log:
        subprocess.run([str(init), "-y", "--no-modify-path", "--profile", "minimal", "--default-host", HOST,
                        "--default-toolchain", state["channel"]], env=env, stdout=log, stderr=subprocess.STDOUT, check=True)
        subprocess.run([str(target / "cargo/bin/rustup"), "component", "add", "rustfmt", "--toolchain", state["channel"]],
                       env=env, stdout=log, stderr=subprocess.STDOUT, check=True)
    path = target / "rustup/toolchains" / (state["channel"] + "-" + HOST) / "bin"
    toolchain_descriptor(path, state)
    atomic_new(target / "owner.json", canonical(owner_value) + b"\n")
    return path


def cold_build(base, source, state, args):
    tools = obtain_toolchain(base, state, args.toolchain_bin, args.offline)
    compiler = toolchain_descriptor(tools, state)
    work = Path(tempfile.mkdtemp(prefix=".build-work-", dir=base))
    for relative, wanted in state["files"].items():
        raw, _ = read_stable(source / relative)
        require(sha(raw) == wanted, "Source changed before cold build")
        target = work / "source" / relative; target.parent.mkdir(parents=True, exist_ok=True); target.write_bytes(raw)
    cargo_home = base / "cargo-cache"; no_links(cargo_home, base); cargo_home.mkdir(exist_ok=True)
    require(cargo_home.is_dir() and not any((cargo_home / name).exists() for name in ("config", "config.toml")),
            "Foreign Cargo cache configuration is refused")
    env = {key: value for key, value in os.environ.items()
           if not key.startswith(("CARGO_", "RUST", "RUSTUP_")) and key not in ("LD_PRELOAD", "LD_LIBRARY_PATH")}
    env.update(CARGO_HOME=str(cargo_home), CARGO_TARGET_DIR=str(work / "target"),
               RUSTC=str(tools / "rustc"), RUSTDOC=str(tools / "rustdoc"),
               RUSTFLAGS="", CARGO_ENCODED_RUSTFLAGS="", CARGO_INCREMENTAL="0")
    cwd = work / "source/rust"
    for parent in (cwd, *cwd.parents):
        require(not any((parent / ".cargo" / name).exists() for name in ("config", "config.toml")),
                "Unpinned ancestor Cargo configuration is refused")
    argv = [str(tools / "cargo"), "build", "--locked", "--bin", "priva", "--jobs", str(args.jobs)]
    if args.build_profile == "release": argv.append("--release")
    if args.offline: argv.append("--offline")
    log_path = work / "build.log"
    with log_path.open("xb") as log:
        subprocess.run(argv, cwd=cwd, env=env, stdout=log, stderr=subprocess.STDOUT, check=True)
    binary = work / "target" / ("release" if args.build_profile == "release" else "debug") / "priva"
    profile = dict(name=args.build_profile, cargo_manifest_sha256=state["files"]["rust/Cargo.toml"],
                   explicit_rustflags="", incremental=False, fast_math=False,
                   cargo_profile_table=state["cargo_profiles"].get(args.build_profile, {}))
    return dict(binary=binary, binary_sha256=digest(binary), evidence={"build.log": read_stable(log_path)[0]},
        profile=profile, toolchain=compiler,
        qualification=dict(origin="cold_locked_build", scientific_acceptance=False, production_backend_selected=False),
        build_command=argv, build_environment={key: env[key] for key in ("CARGO_HOME", "CARGO_TARGET_DIR", "RUSTC", "RUSTDOC", "RUSTFLAGS", "CARGO_ENCODED_RUSTFLAGS", "CARGO_INCREMENTAL")})


def copy_binary(source, target, wanted):
    no_links(source); before = identity(source)
    require(stat.S_ISREG(before["mode"]) and before["mode"] & 0o111, "Expected an executable ordinary source file")
    h = hashlib.sha256()
    with os.fdopen(os.open(source, os.O_RDONLY | os.O_NOFOLLOW), "rb") as src, target.open("xb") as dst:
        require(fd_identity(src.fileno()) == before, "Executable changed before open")
        header = src.read(64)
        require(header[:6] == b"\x7fELF\x02\x01" and header[18:20] == b"\x3e\x00", "Only Linux x86_64 ELF imports are supported")
        dst.write(header); h.update(header)
        for block in iter(lambda: src.read(1024 * 1024), b""):
            dst.write(block); h.update(block)
        dst.flush(); os.fsync(dst.fileno())
    require(identity(source) == before and h.hexdigest() == wanted, "Executable changed or differs from its pinned build receipt")
    target.chmod(0o555)


def install(args):
    require(platform.system() == "Linux" and platform.machine() == "x86_64", "This bootstrap supports Linux x86_64")
    source = Path(args.source_root); state = source_state(source)
    imported = import_evidence(args, state)
    base = managed_base(args.cache_root, create=True)
    with locked(base, True):
        if imported:
            descriptor = artifact_descriptor(state, imported)
            target = base / sha(canonical(descriptor))
            require(args.identity in (None, target.name), "Explicit identity differs from the qualified import")
            if publication_path(target).exists() and not (target / "receipt.json").exists():
                pending, _, _ = load(publication_path(target))
                require(all(pending["members"]["evidence/" + name]["sha256"] == sha(raw)
                            for name, raw in imported["evidence"].items()), "Pending import has different explicit evidence")
                resume_publication(target)
            if target.exists() or target.is_symlink():
                result = verify_artifact(target, state)
                existing, _, _ = load(target / "receipt.json")
                require(existing["evidence_sha256"] == {name: sha(raw) for name, raw in imported["evidence"].items()},
                        "Existing artifact has different import evidence; retained unchanged")
                require(source_state(source) == state, "Source changed during warm import reuse")
                return result
        else:
            ready = find_existing(base, state, args.identity, recover=True)
            if ready:
                require(source_state(source) == state, "Source changed during warm reuse")
                return ready
        spec = imported or cold_build(base, source, state, args)
        descriptor = artifact_descriptor(state, spec); target = base / sha(canonical(descriptor))
        stage = Path(tempfile.mkdtemp(prefix=".artifact-work-", dir=base))
        (stage / "bin").mkdir(); (stage / "evidence").mkdir()
        try:
            copy_binary(spec["binary"], stage / "bin/priva", spec["binary_sha256"])
            version = probe(stage / "bin/priva", state["package_version"])
            if imported:
                require(all(identity(Path(path)) == before for path, before in imported["input_observations"].items()),
                        "Qualified import evidence changed during copy/probe")
            else:
                spec["evidence"]["build.json"] = canonical(dict(artifact=descriptor,
                    command=spec["build_command"], environment=spec["build_environment"])) + b"\n"
            for name, raw in spec["evidence"].items():
                (stage / "evidence" / name).write_bytes(raw)
            save_new(stage / "receipt.json", dict(kind=KIND, schema_version=1, artifact=descriptor,
                evidence_sha256={name: sha(raw) for name, raw in spec["evidence"].items()},
                qualification=spec["qualification"], executable_version=version,
                created_utc=datetime.now(timezone.utc).isoformat(), backend_activated=False))
            require(source_state(source) == state, "Source changed before atomic publication")
            publish_artifact(stage, target)
            result = verify_artifact(target, state)
            require(source_state(source) == state, "Source changed during publication verification")
            return result
        except BaseException as exc:
            if stage.is_dir() and not (stage / "failure.json").exists():
                save_new(stage / "failure.json", dict(error_type=type(exc).__name__, error=str(exc), ready=False))
            raise


def check(args):
    require(platform.system() == "Linux" and platform.machine() == "x86_64", "This bootstrap supports Linux x86_64")
    source = Path(args.source_root); state = source_state(source)
    base = managed_base(args.cache_root)
    with locked(base, False):
        result = find_existing(base, state, args.identity)
        require(result is not None, "No matching managed PriVA Rust CLI")
        require(source_state(source) == state, "Source changed during runtime check")
        return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest="procedure", required=True)
    for name, action in (("install", install), ("check", check)):
        command = sub.add_parser(name)
        command.add_argument("--source-root", type=Path, required=True)
        command.add_argument("--cache-root", type=Path, required=True)
        command.add_argument("--identity")
        command.set_defaults(action=action)
        if name == "install":
            command.add_argument("--qualified-binary", type=Path)
            command.add_argument("--build-pins", type=Path)
            command.add_argument("--build-pins-sha256")
            command.add_argument("--validation", type=Path)
            command.add_argument("--validation-sha256")
            command.add_argument("--toolchain-bin", type=Path)
            command.add_argument("--build-profile", choices=("dev", "release"), default="release")
            command.add_argument("--jobs", type=int, default=1)
            command.add_argument("--offline", action="store_true")
    args = parser.parse_args()
    if args.procedure == "install": require(args.jobs > 0, "--jobs must be positive")
    print(json.dumps(args.action(args), sort_keys=True, allow_nan=False))


if __name__ == "__main__":
    try:
        main()
    except (OSError, ValueError, KeyError, subprocess.SubprocessError) as error:
        sys.exit("PriVA Rust runtime: " + str(error))
