#!/usr/bin/env python3
"""Validate and record official VEP plugin releases without touching other files."""

import argparse
import difflib
import fcntl
import hashlib
import json
import os
from pathlib import Path
import re
import shutil
import sys
import tempfile


REPOSITORY = "https://raw.githubusercontent.com/Ensembl/VEP_plugins"
RECEIPT = ".priva-vep-plugins.json"
UTR_PROFILE = "utrannotator-116-all-effects-20261024-v1"
UTR_UPSTREAM_SHA = "a002fe5a654845627ce2d1f08b4f475db41a051e24002f862f4fc4f4c18d60b0"
UTR_INSTALLED_SHA = "0e4814fb71c9cc1af48a70830ae33434ba620e517f94ed8436979510484e9138"
UTR_REPORTING_MODULE = Path(__file__).resolve().parents[1] / "data/vep/UTRAnnotator116_all_effects.pm"


def digest(path):
    if path.is_symlink() or not path.is_file() or path.stat().st_size == 0:
        raise ValueError(f"Missing or empty regular plugin file: {path}")
    return hashlib.sha256(path.read_bytes()).hexdigest()


def load_receipt(directory):
    try:
        value = json.loads((directory / RECEIPT).read_text())
    except (OSError, ValueError):
        return {}
    return value if isinstance(value, dict) else {}


def local_patch_matches(directory, receipt, filename):
    """Recognize an explicitly recorded patch without calling it upstream code.

    2026-10-24 CT: retain the approved installed UTRAnnotator reporting repair.
    Keep the official hash in plugins; separately check the installed bytes,
    original backup and exact patch. A release/source change needs new review.
    """
    try:
        patch = receipt.get("local_patches", {}).get(filename, {})
        return bool(
            patch.get("upstream_sha256") == receipt["plugins"][filename]
            and re.fullmatch(r"[0-9a-f]{64}", patch.get("installed_sha256", ""))
            and re.fullmatch(r"[0-9a-f]{64}", patch.get("patch_sha256", ""))
            and patch["installed_sha256"] == digest(directory / filename)
            and patch["upstream_sha256"] == digest(Path(patch["original_file"]))
            and patch["patch_sha256"] == digest(Path(patch["patch_file"]))
            and bool(patch.get("profile"))
        )
    except (OSError, ValueError, TypeError, KeyError, AttributeError):
        return False


def matches(directory, receipt, release, names):
    if (receipt.get("schema_version") != 1 or receipt.get("release") != release
            or receipt.get("source") != f"{REPOSITORY}/release/{release}"
            or not re.fullmatch(r"[0-9a-f]{64}", str(receipt.get("plugin_config_sha256", "")))
            or not isinstance(receipt.get("plugins"), dict)):
        return False
    return all(
        local_patch_matches(directory, receipt, name + ".pm")
        if name + ".pm" in receipt.get("local_patches", {}) else
        receipt["plugins"].get(name + ".pm") == digest(directory / (name + ".pm"))
        for name in names)


def write_receipt(path, value):
    with tempfile.NamedTemporaryFile(mode="w", dir=path.parent, prefix=".priva-plugin-receipt-",
                                     delete=False, encoding="utf-8") as handle:
        temporary = Path(handle.name)
        json.dump(value, handle, indent=2, sort_keys=True)
        handle.write("\n")
    try:
        os.replace(temporary, path)
    finally:
        temporary.unlink(missing_ok=True)


def ensure_utr_reporting(directory, release):
    """Reproduce only the approved VEP116 reporting profile, without recalculation.

    User, 2026-10-24 CT: "So maybe you can modify the pm file directly installed
    in this cluster ?" The approved change reports all already-calculated UTR
    effects deterministically. It does not change the five biological methods.
    2026-10-26 CT: make those exact installed bytes reproducible by the installer;
    never infer compatibility with another upstream source or release. Existing
    recorded installations keep their original backup paths and receipt bytes.
    """
    if release != "116":
        return {"status": "not_applicable", "release": release}
    target, receipt_path = directory / "UTRAnnotator.pm", directory / RECEIPT
    pending = directory / ".priva-utr-reporting.pending.json"
    lock = directory / ".priva-utr-reporting.lock"
    if any(path.is_symlink() for path in (target, receipt_path, pending, lock)):
        raise ValueError("Refusing symlinked UTR reporting installation files")
    # A valid completed installation is a read-only no-op, including on this
    # cluster's existing recorded patch. Do not create a new lock/backup there.
    prior = load_receipt(directory)
    patches = prior.get("local_patches", {})
    if (prior.get("release") == release
            and prior.get("source") == f"{REPOSITORY}/release/{release}"
            and isinstance(prior.get("plugins"), dict)
            and prior["plugins"].get("UTRAnnotator.pm") == UTR_UPSTREAM_SHA
            and isinstance(patches, dict)
            and isinstance(patches.get("UTRAnnotator.pm"), dict)
            and patches["UTRAnnotator.pm"].get("profile") == UTR_PROFILE
            and patches["UTRAnnotator.pm"].get("installed_sha256") == UTR_INSTALLED_SHA
            and not pending.exists()
            and local_patch_matches(directory, prior, "UTRAnnotator.pm")):
        return {"status": "reused", "profile": UTR_PROFILE}
    with lock.open("a") as guard:
        fcntl.flock(guard, fcntl.LOCK_EX)
        receipt = load_receipt(directory)
        if (receipt.get("release") != release
                or receipt.get("source") != f"{REPOSITORY}/release/{release}"
                or not isinstance(receipt.get("plugins"), dict)
                or receipt["plugins"].get("UTRAnnotator.pm") != UTR_UPSTREAM_SHA):
            raise ValueError("UTR reporting fix requires the pinned official VEP116 module; files were preserved")
        patches = receipt.get("local_patches", {})
        if not isinstance(patches, dict) or ("UTRAnnotator.pm" in patches
                                            and not isinstance(patches["UTRAnnotator.pm"], dict)):
            raise ValueError("Malformed UTR reporting patch receipt; files were preserved")
        recorded = patches.get("UTRAnnotator.pm")
        if recorded:
            if (recorded.get("profile") != UTR_PROFILE
                    or recorded.get("installed_sha256") != UTR_INSTALLED_SHA
                    or not local_patch_matches(directory, receipt, "UTRAnnotator.pm")):
                raise ValueError("Existing UTR reporting patch differs from the approved profile; files were preserved")
            if pending.exists():
                recovery = json.loads(pending.read_text())
                if (not isinstance(recovery, dict)
                        or set(recovery) != {"receipt_before", "receipt_after"}
                        or recovery["receipt_after"] != receipt
                        or not isinstance(recovery["receipt_before"], dict)):
                    raise ValueError("Conflicting UTR reporting recovery receipt; files were preserved")
                before = recovery["receipt_before"]
                old_patches = before.get("local_patches", {})
                if not isinstance(old_patches, dict) or "UTRAnnotator.pm" in old_patches:
                    raise ValueError("Conflicting UTR reporting recovery origin; files were preserved")
                expected = json.loads(json.dumps(before))
                expected.setdefault("local_patches", {})["UTRAnnotator.pm"] = recorded
                if expected != receipt:
                    raise ValueError("Conflicting UTR reporting recovery origin; files were preserved")
                pending.unlink()
            return {"status": "reused", "profile": UTR_PROFILE}
        before_sha = digest(receipt_path)
        observed_sha = digest(target)
        if observed_sha not in (UTR_UPSTREAM_SHA, UTR_INSTALLED_SHA):
            raise ValueError("UTRAnnotator differs from the pinned original; files were preserved")
        recovery = json.loads(pending.read_text()) if pending.exists() else None
        if pending.exists() and (not isinstance(recovery, dict)
                                 or set(recovery) != {"receipt_before", "receipt_after"}):
            raise ValueError("Malformed UTR reporting recovery receipt; files were preserved")
        if observed_sha == UTR_INSTALLED_SHA and recovery is None:
            raise ValueError("Unrecorded UTR reporting customization; files were preserved")
        # This is bundled executable source, not a downloaded annotation table.
        # Its output hash proves identity to the already approved plugin profile.
        candidate = UTR_REPORTING_MODULE.read_bytes()
        if hashlib.sha256(candidate).hexdigest() != UTR_INSTALLED_SHA:
            raise ValueError("Bundled UTR reporting implementation does not reproduce the approved profile")
        folder = directory / ".priva-local-patches" / "UTRAnnotator116_all_effects"
        if folder.is_symlink() or folder.parent.is_symlink():
            raise ValueError("Refusing a symlinked UTR patch backup directory")
        original, patch_file = folder / "UTRAnnotator.original.pm", folder / "UTRAnnotator_all_effects.patch"
        upstream = target.read_bytes() if observed_sha == UTR_UPSTREAM_SHA else original.read_bytes()
        if hashlib.sha256(upstream).hexdigest() != UTR_UPSTREAM_SHA:
            raise ValueError("UTR reporting recovery original differs; files were preserved")
        patch = "".join(difflib.unified_diff(
            upstream.decode().splitlines(True), candidate.decode().splitlines(True),
            fromfile="UTRAnnotator.original.pm", tofile="UTRAnnotator.pm")).encode()
        metadata = {"profile": UTR_PROFILE, "upstream_sha256": UTR_UPSTREAM_SHA,
                    "installed_sha256": UTR_INSTALLED_SHA,
                    "original_file": str(original), "patch_file": str(patch_file),
                    "patch_sha256": hashlib.sha256(patch).hexdigest()}
        updated = json.loads(json.dumps(receipt))
        updated.setdefault("local_patches", {})["UTRAnnotator.pm"] = metadata
        transaction = {"receipt_before": receipt, "receipt_after": updated}
        if recovery is not None and recovery != transaction:
            raise ValueError("UTR reporting recovery context changed; files were preserved")
        folder.mkdir(parents=True, exist_ok=True)
        for path, contents in ((original, upstream), (patch_file, patch)):
            if path.is_symlink() or (path.exists() and path.read_bytes() != contents):
                raise ValueError(f"Conflicting UTR reporting backup; file was preserved: {path}")
            if not path.exists():
                with path.open("xb") as handle:
                    handle.write(contents)
        # Record the exact recoverable transaction before either live write.
        write_receipt(pending, transaction)
        if digest(receipt_path) != before_sha or digest(target) != observed_sha:
            raise ValueError("UTR reporting files changed during preparation; live files were preserved")
        if observed_sha == UTR_UPSTREAM_SHA:
            with tempfile.NamedTemporaryFile(dir=directory, prefix=".priva-utr-", delete=False) as handle:
                temporary = Path(handle.name)
                handle.write(candidate)
            try:
                temporary.chmod(target.stat().st_mode & 0o777)
                os.replace(temporary, target)
            finally:
                temporary.unlink(missing_ok=True)
        write_receipt(receipt_path, updated)
        if not local_patch_matches(directory, load_receipt(directory), "UTRAnnotator.pm"):
            raise ValueError("UTR reporting publication is incomplete; recovery receipt was retained")
        pending.unlink()
        return {"status": "installed", "profile": UTR_PROFILE}


def capture_transition(path, config, previous, release):
    if not re.fullmatch(r"[1-9][0-9]*", previous or "") or previous == release:
        raise ValueError("A release transition needs distinct numeric source and target releases")
    if not isinstance(config, dict):
        raise ValueError("The installer configuration must be a mapping")
    root = config.get("cache_root")
    if root is not None and root != "":
        if not isinstance(root, str) or not re.fullmatch(r"/[A-Za-z0-9_./-]+", root):
            raise ValueError("cache_root must be an absolute path without spaces or shell characters")
        root = os.path.normpath(root)
        if not root.strip("/"):
            raise ValueError("cache_root must name a dedicated cache directory")
        folders = config.get("cache_folders") or {}
        if not isinstance(folders, dict):
            raise ValueError("cache_folders must be a mapping")
        folder = folders.get("vep_plugin", "VEP_plugin")
        if not isinstance(folder, str) or not re.fullmatch(r"[A-Za-z0-9_.-]+", folder) or folder in (".", ".."):
            raise ValueError("cache_folders.vep_plugin must be a single folder name")
        directory = Path(root) / folder
    else:
        directory = Path(config.get("vep_plugins_dir") or "")
    if not directory.is_absolute():
        raise ValueError("vep_plugins_dir must resolve to an absolute path before a release switch")
    directory = directory.resolve()
    prior = {}
    for source in directory.glob("*.pm"):
        if re.fullmatch(r"[A-Za-z][A-Za-z0-9_]*\.pm", source.name) and not source.is_symlink():
            prior[source.name] = digest(source)
    write_receipt(path, {
        "schema_version": 1, "from_release": previous, "to_release": release,
        "environment_prefix": str(path.parent.resolve()),
        "plugins_directory": str(directory), "plugins": prior,
    })


def matching_transition(path, directory, release):
    if path is None:
        return {}
    try:
        value = json.loads(path.read_text())
        if (not isinstance(value, dict) or value.get("schema_version") != 1
                or value.get("to_release") != release or value.get("from_release") == release
                or not re.fullmatch(r"[1-9][0-9]*", str(value.get("from_release", "")))
                or value.get("environment_prefix") != str(path.parent.resolve())
                or value.get("plugins_directory") != str(directory.resolve())
                or not isinstance(value.get("plugins"), dict)):
            return {}
        return value
    except (OSError, ValueError):
        return {}


def install(directory, release, names, stage, transition_file=None):
    allow_release_change = False
    transition = matching_transition(transition_file, directory, release)
    # Validate the whole download before considering any shared-file writes.
    config_sha256 = digest(stage / "plugin_config.txt")
    hashes = {name + ".pm": digest(stage / (name + ".pm")) for name in names}
    for name in names:
        if not re.search(r"(?m)^\s*package\s+" + re.escape(name) + r"\s*;",
                         (stage / (name + ".pm")).read_text()):
            raise ValueError(f"Downloaded plugin does not declare package {name}: {stage / (name + '.pm')}")
    previous = load_receipt(directory)
    old_release = previous.get("release", "")
    old_plugins = previous.get("plugins", {})
    local_patches = previous.get("local_patches", {})
    if not isinstance(local_patches, dict):
        raise ValueError("Malformed local plugin patch receipt; installed files were preserved")
    for filename, patch in local_patches.items():
        if (not local_patch_matches(directory, previous, filename)
                or old_release != release
                or previous.get("source") != f"{REPOSITORY}/release/{release}"
                or filename not in hashes or hashes[filename] != patch["upstream_sha256"]):
            raise ValueError(
                f"Recorded local plugin patch needs review before an upstream change: {filename}. "
                "Installed files were preserved; restore the recorded original before changing releases."
            )
    if (isinstance(old_release, str) and re.fullmatch(r"[1-9][0-9]*", old_release)
            and old_release != release and isinstance(old_plugins, dict)):
        try:
            old_names = [name[:-3] for name in old_plugins if re.fullmatch(r"[A-Za-z][A-Za-z0-9_]*\.pm", name)]
            # A failed publication can leave selected files at the new release.
            # Retry only if old files match their recorded hash or this download.
            metadata_valid = (previous.get("schema_version") == 1
                              and previous.get("source") == f"{REPOSITORY}/release/{old_release}"
                              and re.fullmatch(r"[0-9a-f]{64}", str(previous.get("plugin_config_sha256", ""))))
            known_files = all(
                digest(directory / (name + ".pm")) in {old_plugins[name + ".pm"], hashes.get(name + ".pm")}
                for name in old_names
            )
            allow_release_change |= bool(metadata_valid and old_names and known_files)
        except (OSError, ValueError):
            pass
    for filename, expected in hashes.items():
        destination = directory / filename
        if destination.is_symlink():
            raise ValueError(f"Refusing to replace a symlinked VEP plugin: {destination}")
        if (filename not in local_patches and destination.exists() and digest(destination) != expected and not allow_release_change
                and transition.get("plugins", {}).get(filename) != digest(destination)):
            raise ValueError(
                f"Existing VEP plugin differs from official release {release}: {destination}. "
                "Choose an empty vep_plugins_dir for this release or explicitly switch the environment's major release; "
                "existing customized or unrecorded plugins were preserved."
            )
    receipt = {
        "schema_version": 1,
        "release": release,
        "source": f"{REPOSITORY}/release/{release}",
        "plugin_config_sha256": config_sha256,
        "plugins": hashes,
    }
    if local_patches:
        receipt["local_patches"] = local_patches
    for filename, expected in hashes.items():
        destination = directory / filename
        if filename in local_patches or (destination.exists() and digest(destination) == expected):
            continue
        with tempfile.NamedTemporaryFile(dir=directory, prefix=".priva-plugin-", delete=False) as handle:
            temporary = Path(handle.name)
        try:
            shutil.copy2(stage / filename, temporary)
            os.replace(temporary, destination)
        finally:
            temporary.unlink(missing_ok=True)
    write_receipt(directory / RECEIPT, receipt)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("mode", choices=("check", "install", "transition", "clear-transition", "ensure-utr-reporting"))
    parser.add_argument("--directory", type=Path)
    parser.add_argument("--release", required=True)
    parser.add_argument("--plugins", default="")
    parser.add_argument("--stage", type=Path)
    parser.add_argument("--transition-file", type=Path)
    parser.add_argument("--from-release")
    parser.add_argument("--config-json")
    args = parser.parse_args()
    try:
        if not re.fullmatch(r"[1-9][0-9]*", args.release):
            raise ValueError("VEP release must be a numeric major release")
        if args.mode == "transition":
            if args.transition_file is None:
                raise ValueError("A release transition requires an environment transition file")
            capture_transition(args.transition_file, json.loads(args.config_json or "null"), args.from_release, args.release)
            return 0
        if args.directory is None:
            raise ValueError("A plugin directory is required")
        if args.mode == "clear-transition":
            if matching_transition(args.transition_file, args.directory, args.release):
                args.transition_file.unlink()
            return 0
        if args.mode == "ensure-utr-reporting":
            if not args.directory.is_dir():
                raise ValueError(f"VEP plugin directory does not exist: {args.directory}")
            print(json.dumps(ensure_utr_reporting(args.directory, args.release), sort_keys=True))
            return 0
        names = list(dict.fromkeys(args.plugins.split(",")))
        if not names or any(not re.fullmatch(r"[A-Za-z][A-Za-z0-9_]*", name) for name in names):
            raise ValueError("VEP plugin names must be simple module names")
        if not args.directory.is_dir():
            raise ValueError(f"VEP plugin directory does not exist: {args.directory}")
        if args.mode == "check":
            if not matches(args.directory, load_receipt(args.directory), args.release, names):
                return 1
        else:
            if args.stage is None:
                raise ValueError("Installing plugin files requires a completed staging directory")
            install(args.directory, args.release, names, args.stage, args.transition_file)
    except (OSError, ValueError, TypeError) as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
