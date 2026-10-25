#!/usr/bin/env python3
"""Publish a generated motif/PM1 pickle and its existing mapping receipt together.

``run`` changes only the scientific command's --output argument. Failed work is
retained in a private same-directory stage. Only a completed staged pair may
create a pending marker; ``recover`` can finish its two recorded renames before
the installer's ordinary artifact-route checks. An unknown unreceipted artifact
is never recovery evidence. This is process-interruption recovery, not a new
scientific cache-validity or bundled-snapshot policy.
"""
from __future__ import annotations

import argparse
from contextlib import contextmanager
import fcntl
import importlib.util
import json
import os
from pathlib import Path
import stat
import subprocess
import sys
import tempfile


SPEC = importlib.util.spec_from_file_location(
    "mapped_artifact_biomart", Path(__file__).with_name("biomart_nm_to_enst.py"))
BIOMART = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(BIOMART)
PRODUCERS = {"motifs": "am_pick_missense_constrained_segments.py", "pm1": "prepare_pm1_regions.py"}
SCHEMA = "priva-mapped-artifact-publication-v1"


def require(condition, message):
    if not condition:
        raise ValueError(message)


def receipt_path(artifact):
    return Path(str(artifact) + ".nm_to_enst.json")


def pending_path(artifact):
    return Path(str(artifact) + ".publish.pending.json")


def file_state(path):
    """Identity plus content for these generated artifacts, never a cache scan."""
    path = Path(path)
    if not os.path.lexists(path):
        return None
    before = path.lstat()
    require(stat.S_ISREG(before.st_mode) and before.st_nlink == 1,
            "Publication files must be regular, unlinked files: " + str(path))
    digest = BIOMART.sha256(path)
    after = path.lstat()
    fields = ("st_ino", "st_size", "st_mtime_ns", "st_ctime_ns", "st_mode", "st_nlink", "st_uid")
    require(all(getattr(before, key) == getattr(after, key) for key in fields),
            "Publication file changed while reading: " + str(path))
    return {"inode": before.st_ino, "size": before.st_size, "mtime_ns": before.st_mtime_ns,
            "ctime_ns": before.st_ctime_ns, "mode": before.st_mode,
            "nlink": before.st_nlink, "uid": before.st_uid, "sha256": digest}


def moved_state_matches(actual, expected):
    # Rename may change ctime; the exact staged inode/content/mode must survive.
    return actual is not None and expected is not None and all(
        actual[key] == value for key, value in expected.items() if key != "ctime_ns")


def context(args):
    producer = Path(args.producer)
    require(producer.is_absolute() and producer.name == PRODUCERS[args.kind],
            "Use the selected motif/PM1 scientific producer")
    require(not producer.is_symlink() and producer.is_file(), "Scientific producer is not a regular file")
    config = Path(args.config)
    require(config.is_absolute() and not config.is_symlink() and config.is_file(),
            "Use the actual regular runtime config")
    result = {"kind": args.kind, "assembly": args.assembly, "uid": os.getuid(),
            "cache_root": str(args.cache_root.resolve()), "output_dir": str(args.output_dir.resolve()),
            "producer": str(producer),
            "producer_sha256": BIOMART.sha256(producer), "config": str(config),
            "config_sha256": BIOMART.sha256(config),
            "helper_sha256": BIOMART.sha256(Path(__file__)),
            "mapping_helper_sha256": BIOMART.sha256(Path(BIOMART.__file__)),
            "mapping_identity": BIOMART.mapping_identity(args.mapping, args.input)}
    if getattr(args, "annotation_base", None) is not None:
        require(args.kind == "motifs", "Only motifs consume the complete AlphaMissense base")
        spec = importlib.util.spec_from_file_location("mapped_alpha_base", producer.with_name("alphamissense_annotation_base.py"))
        base = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(base)
        result["scientific_input"] = base.motif_binding(args.config, producer.parents[1], args.annotation_base)
    if getattr(args, "native_lineage", None) is not None:
        result["native_input_lineage"] = {"path": str(args.native_lineage.absolute()),
                                          "sha256": BIOMART.sha256(args.native_lineage)}
    return result


def artifact_record(expected, artifact_state):
    result = {"mapping_identity": expected["mapping_identity"],
              "artifact_bytes": artifact_state["size"], "artifact_sha256": artifact_state["sha256"]}
    if "scientific_input" in expected:
        result["scientific_input"] = expected["scientific_input"]
    if "native_input_lineage" in expected:
        result["native_input_lineage"] = expected["native_input_lineage"]
    return result


def recovery_context(args, target):
    """Normal installers can finish a pinned, explicit native adoption too."""
    expected = context(args)
    value = load_pending(pending_path(target))
    pinned = value.get("context", {}).get("native_input_lineage")
    if pinned is not None and "native_input_lineage" not in expected:
        require(args.kind == "motifs" and "scientific_input" in expected and
                isinstance(pinned, dict) and set(pinned) == {"path", "sha256"},
                "Unexpected native adoption context")
        path = Path(pinned["path"])
        require(path.is_absolute() and not path.is_symlink() and path.is_file() and
                BIOMART.sha256(path) == pinned["sha256"], "Pinned native lineage changed")
        lineage = json.loads(path.read_text())
        require(lineage.get("schema") == "priva-alphamissense-motif-input-v1" and lineage.get("provenance") and
                lineage.get("scientific_input") == expected["scientific_input"] and
                lineage.get("mapping_identity") == expected["mapping_identity"] and
                lineage.get("artifact_sha256") == value["new"]["artifact"]["sha256"] and
                lineage.get("artifact_bytes") == value["new"]["artifact"]["size"],
                "Pinned native lineage no longer describes the pending motif")
        original = Path(lineage["artifact_path"])
        BIOMART.check_receipt(original, args.mapping, args.input)
        require(file_state(original)["sha256"] == lineage["artifact_sha256"] and
                file_state(receipt_path(original))["sha256"] == lineage["receipt_sha256"],
                "Preserved historical native pair changed")
        expected["native_input_lineage"] = pinned
    return expected


def managed_target(args):
    artifact = Path(args.artifact)
    require(artifact.is_absolute() and ".." not in artifact.parts,
            "Use an absolute selected artifact path")
    return BIOMART._managed_artifact_path(
        cache_root=args.cache_root, output_dir=artifact.parent, relative=Path(artifact.name))


@contextmanager
def directory_lock(directory):
    descriptor = os.open(directory, os.O_RDONLY | os.O_DIRECTORY)
    try:
        fcntl.flock(descriptor, fcntl.LOCK_EX | fcntl.LOCK_NB)
        yield
    finally:
        os.close(descriptor)


def load_pending(marker):
    before = file_state(marker)
    require(before is not None and before["size"] <= 65536, "Missing/oversized pending publication")
    value = json.loads(marker.read_text())
    require(file_state(marker) == before, "Pending publication changed while reading")
    require(isinstance(value, dict) and set(value) ==
            {"schema", "context", "target", "stage", "prior", "new", "producer_job_id"} and
            value["schema"] == SCHEMA, "Unknown pending publication format")
    return value


def finish_pending(args, target, expected_context):
    """Accept only the three states produced by the ordered two-file commit."""
    marker = pending_path(target)
    value = load_pending(marker)
    require(value["context"] == expected_context and value["target"] == str(target),
            "Pending publication producer/config/mapping/target differs")
    stage = Path(value["stage"])
    require(stage.parent == target.parent and stage.name.startswith(".mapped-artifact-" + args.kind + "-") and
            not stage.is_symlink() and stage.is_dir() and stage.stat().st_uid == os.getuid() and
            # A shared scratch parent can legitimately contribute its setgid bit.
            stat.S_IMODE(stage.stat().st_mode) & 0o777 == 0o700,
            "Pending publication stage is not owned beside target")
    require(set(value["prior"]) == {"artifact", "receipt"} and
            set(value["new"]) == {"artifact", "receipt"}, "Pending pair fields differ")
    paths = {"artifact": (stage / target.name, target),
             "receipt": (receipt_path(stage / target.name), receipt_path(target))}
    locations = []
    for key, (staged, final) in paths.items():
        staged_state, final_state = file_state(staged), file_state(final)
        expected = value["new"][key]
        require(isinstance(expected, dict) and set(expected) ==
                {"inode", "size", "mtime_ns", "ctime_ns", "mode", "nlink", "uid", "sha256"} and
                expected.get("size", 0) > 0,
                "Pending staged file identity is incomplete")
        if staged_state is not None:
            require(moved_state_matches(staged_state, expected) and final_state == value["prior"][key],
                    "Pending stage or previous selected artifact changed")
            locations.append("stage")
        else:
            require(moved_state_matches(final_state, expected), "Pending published file changed or disappeared")
            locations.append("final")
    require(tuple(locations) in (("stage", "stage"), ("final", "stage"), ("final", "final")),
            "Pending publication is not an ordered two-file state")
    # The receipt's original contract remains authoritative in a mixed state too.
    staged_receipt, final_receipt = paths["receipt"]
    record = json.loads((staged_receipt if locations[1] == "stage" else final_receipt).read_text())
    require(record == artifact_record(expected_context, value["new"]["artifact"]),
            "Pending mapping receipt does not describe the staged artifact")
    for location, (staged, final) in zip(locations, paths.values()):
        if location == "stage":
            os.replace(staged, final)
    BIOMART.check_receipt(target, args.mapping, args.input)
    require(all(moved_state_matches(file_state(paths[key][1]), value["new"][key]) for key in paths),
            "Published pair changed before commit")
    marker.unlink()
    # The motif CLI's ordinary log stays in this retained stage if present.
    if not any(stage.iterdir()):
        stage.rmdir()


def recover(args):
    # A changed mapping can route the completed pair before YAML names it.
    target = Path(args.artifact)
    if not os.path.lexists(pending_path(target)):
        if not (os.path.lexists(target) or os.path.lexists(receipt_path(target))):
            return False
        if not args.cache_root.is_absolute() or not args.output_dir.is_absolute():
            return False  # An ordinary legacy/supplied warm read is not a publication.
        identity = BIOMART.mapping_identity(args.mapping, args.input)
        candidate = BIOMART.changed_mapping_target(
            identity=identity, cache_root=args.cache_root, output_dir=args.output_dir,
            assembly=args.assembly, kind=args.kind)
        if not os.path.lexists(pending_path(candidate)):
            return False
        require(BIOMART._artifact_receipt(target)["mapping_identity"] != identity,
                "Pending mapping namespace is not the selected artifact's changed mapping")
        target = candidate
    else:
        target = managed_target(args)
    with directory_lock(target.parent):
        if os.path.lexists(pending_path(target)):
            finish_pending(args, target, recovery_context(args, target))
            print("Recovered completed mapped artifact publication: " + str(args.artifact), flush=True)
            return True
    return False


def run(args):
    """Run the unchanged producer against a private output, then commit its pair."""
    target = managed_target(args)
    target.parent.mkdir(parents=True, exist_ok=True)
    with directory_lock(target.parent):
        expected = context(args)
        require(expected["mapping_identity"] == json.loads(args.expected_identity),
                "Selected mapping changed before generation")
        if os.path.lexists(pending_path(target)):
            finish_pending(args, target, recovery_context(args, target))
            return
        if os.path.lexists(target) or os.path.lexists(receipt_path(target)):
            BIOMART.check_receipt(target, args.mapping, args.input)
        prior = {"artifact": file_state(target), "receipt": file_state(receipt_path(target))}
        argv = list(args.arguments)
        if argv and argv[0] == "--":
            argv.pop(0)
        require(argv.count("--output") == 1, "Scientific argv must contain exactly one --output")
        index = argv.index("--output") + 1
        require(index < len(argv) and argv[index] == str(args.artifact),
                "Scientific --output must be the selected configured artifact")
        stage = Path(tempfile.mkdtemp(prefix=".mapped-artifact-" + args.kind + "-", dir=target.parent))
        staged = stage / target.name
        argv[index] = str(staged)
        print("Generating mapped artifact in retained stage: " + str(stage), flush=True)
        subprocess.run([sys.executable, str(args.producer), *argv], check=True)
        staged_state = file_state(staged)
        require(staged_state is not None and staged_state["size"] > 0,
                "Scientific producer created no nonempty pickle")
        require(context(args) == expected, "Producer/config/mapping changed during generation")
        BIOMART.atomic_json(receipt_path(staged), artifact_record(expected, staged_state))
        BIOMART.check_receipt(staged, args.mapping, args.input)
        require({"artifact": file_state(target), "receipt": file_state(receipt_path(target))} == prior,
                "Selected artifact changed during generation")
        value = {"schema": SCHEMA, "context": expected, "target": str(target), "stage": str(stage),
                 "prior": prior, "new": {"artifact": file_state(staged),
                                           "receipt": file_state(receipt_path(staged))},
                 "producer_job_id": os.environ.get("SLURM_JOB_ID")}
        temporary = stage / "publication.json"
        with temporary.open("x") as handle:
            json.dump(value, handle, indent=2, sort_keys=True)
            handle.write("\n")
        require(not os.path.lexists(pending_path(target)), "Another pending publication exists")
        os.replace(temporary, pending_path(target))
        finish_pending(args, target, expected)


def adopt_native(args, source, lineage):
    """Publish a verified copy with a new bound receipt; retain history intact."""
    source = Path(source)
    target = managed_target(args)
    require(source.resolve() != target.resolve(), "Native motif adoption must select a new managed target")
    require(not source.is_symlink() and source.is_file(), "Preserved native motif must be an independent regular file")
    target.parent.mkdir(parents=True, exist_ok=True)
    with directory_lock(target.parent):
        expected = context(args)
        require(lineage.get("scientific_input") == expected.get("scientific_input"),
                "Preserved native motif base/producer binding differs")
        if os.path.lexists(pending_path(target)):
            finish_pending(args, target, expected)
            return
        if os.path.lexists(target) or os.path.lexists(receipt_path(target)):
            BIOMART.check_receipt(target, args.mapping, args.input)
            completed = json.loads(receipt_path(target).read_text())
            require(completed.get("scientific_input") == expected.get("scientific_input") and
                    completed.get("native_input_lineage") == expected.get("native_input_lineage") and
                    completed.get("artifact_sha256") == lineage.get("artifact_sha256"),
                    "Native adoption must not replace an existing artifact or historical receipt")
            return
        BIOMART.check_receipt(source, args.mapping, args.input)
        source_before, receipt_before = file_state(source), file_state(receipt_path(source))
        original = json.loads(receipt_path(source).read_text())
        require(all(lineage.get(key) == original.get(key) for key in
                    ("mapping_identity", "artifact_bytes", "artifact_sha256")), "Preserved native motif identity differs")
        require(lineage.get("receipt_sha256") == receipt_before["sha256"], "Preserved native receipt SHA256 differs")
        stage = Path(tempfile.mkdtemp(prefix=".mapped-artifact-motifs-", dir=target.parent))
        staged = stage / target.name
        # No hardlinks: the historical first result remains physically intact.
        with source.open("rb") as incoming, staged.open("xb") as output:
            import shutil
            shutil.copyfileobj(incoming, output, 8 * 1024 * 1024)
        staged_state = file_state(staged)
        require(staged_state["sha256"] == source_before["sha256"] and
                (file_state(source), file_state(receipt_path(source))) == (source_before, receipt_before),
                "Preserved native motif changed while copying")
        BIOMART.atomic_json(receipt_path(staged), artifact_record(expected, staged_state))
        require(context(args) == expected, "Native motif inputs changed before publication")
        value = {"schema": SCHEMA, "context": expected, "target": str(target), "stage": str(stage),
                 "prior": {"artifact": None, "receipt": None},
                 "new": {"artifact": file_state(staged), "receipt": file_state(receipt_path(staged))},
                 "producer_job_id": os.environ.get("SLURM_JOB_ID")}
        BIOMART.atomic_json(pending_path(target), value)
        finish_pending(args, target, expected)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    commands = parser.add_subparsers(dest="command", required=True)
    for name, function in (("recover", recover), ("run", run)):
        child = commands.add_parser(name)
        child.add_argument("--artifact", type=Path, required=True)
        child.add_argument("--mapping", type=Path, required=True)
        child.add_argument("--input", type=Path, action="append", default=[])
        child.add_argument("--cache-root", type=Path, required=True)
        child.add_argument("--output-dir", type=Path, required=True)
        child.add_argument("--assembly", required=True)
        child.add_argument("--config", type=Path, required=True)
        child.add_argument("--producer", type=Path, required=True)
        child.add_argument("--kind", choices=PRODUCERS, required=True)
        child.add_argument("--annotation-base", type=Path,
                           help="Bind motif publication atomically to its complete annotation and scientific producer")
        child.add_argument("--native-lineage", type=Path,
                           help="Exact closed native lineage when recovering an explicitly adopted motif pair")
        if name == "run":
            child.add_argument("--expected-identity", required=True)
            child.add_argument("arguments", nargs=argparse.REMAINDER)
        child.set_defaults(function=function)
    args = parser.parse_args()
    try:
        args.function(args)
    except (OSError, ValueError, subprocess.CalledProcessError) as error:
        print("ERROR: mapped artifact publication: " + str(error), file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
