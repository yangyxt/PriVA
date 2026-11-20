#!/usr/bin/env python3
"""Compare `priva vep-plugins` with its unchanged Python source.

check --binary PRIVA --repo CHECKOUT --work-root NEW_DIRECTORY [--only CASE ...]
oracle --producer SCRIPT --exception JSON -- ARGUMENTS (internal reference launcher)

Each case builds two identical fixture folders (plugins, stage, environment
prefix). Setup steps that need an earlier installation run the unchanged
script itself in each folder before the comparison, so recorded patches,
backups and recovery receipts are the ones it writes. The Python arm then runs
CHECKOUT/scripts/vep_plugins_release.py through the oracle launcher, which
records the class and text of the exception main() reports without changing
main() or any check; the native arm runs `PRIVA vep-plugins` with the same
arguments and PRIVA_DATA_DIR=CHECKOUT/data. The UTRAnnotator 116 cases use
the pinned upstream module in tests/fixtures.

Compared exactly: exit status, stdout, the exception category and message
(case folder written as <ARM>, temporary-file suffixes as <TEMP>), and the
complete fixture tree after the run: every entry's kind, permission bits and
content (text with the case folder written as <ARM>), and for each plugin
file whether its modification time equals the staged file's. Argument errors
exit 2 in both; argparse's usage text is not compared. No case reaches the
network. Writes WORK_ROOT/complete.json and exits 1 on any mismatch.
"""
from __future__ import annotations

import argparse
import hashlib
import importlib.util
import json
import os
from pathlib import Path
import re
import shutil
import subprocess
import sys

sys.dont_write_bytecode = True
PRODUCER_SHA256 = "a16f8a49aded579807f0e1758bac859028452e80cd5d1ac5de3e3956e5e09254"
NAMES = "SpliceAI,SpliceVault"
UTR = "UTRAnnotator"
STAGED_NS = 1_000_000_000_123_456_789
TEMPORARY = re.compile(r"(\.priva-(?:plugin-receipt-|plugin-|utr-))[A-Za-z0-9_]+")


def load(path, name):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def write_json(path, value):
    Path(path).write_text(json.dumps(value, sort_keys=True, indent=2) + "\n")


def oracle(args):
    """Run the unchanged main() and record the exception it reports."""
    module = load(args.producer, "original_vep_plugins_release")
    record = {}
    builtin_print = print

    def observed_print(*values, **options):
        # main() reports a caught failure with print(f"ERROR: {exc}") inside
        # its except clause, where sys.exc_info() is that exception.
        if options.get("file") is sys.stderr and values and str(values[0]).startswith("ERROR: "):
            kind, value, _ = sys.exc_info()
            if kind is not None:
                record.update(category=kind.__name__, message=str(value))
        return builtin_print(*values, **options)

    module.print = observed_print
    sys.argv = [str(args.producer), *args.arguments]
    try:
        status = module.main()
    except SystemExit:
        raise
    except BaseException as exc:
        record.update(category=type(exc).__name__, message=str(exc), uncaught=True)
        write_json(args.exception, record)
        raise
    write_json(args.exception, record)
    raise SystemExit(status)


# ---------------------------------------------------------------- fixtures

def module(name, extra=""):
    return f"package {name};\n# fixture{extra}\n1;\n"


def stage_files(release, names=NAMES, utr=False):
    files = {f"stage{release}/plugin_config.txt": "official config fixture " + release}
    for name in names.split(","):
        files[f"stage{release}/{name}.pm"] = module(name, " " + release)
    if utr:
        files[f"stage{release}/{UTR}.pm"] = ("copy", "tests/fixtures/UTRAnnotator116.original.pm")
    return files


def materialize(arm, case, repo):
    arm.mkdir(parents=True)
    for folder in ("plugins", "env"):
        (arm / folder).mkdir()
    for relative, spec in case.get("files", {}).items():
        path = arm / relative
        path.parent.mkdir(parents=True, exist_ok=True)
        if isinstance(spec, str):
            path.write_text(spec.replace("{arm}", str(arm)))
        elif spec[0] == "bytes":
            path.write_bytes(spec[1])
        elif spec[0] == "copy":
            shutil.copyfile(repo / spec[1], path)
        elif spec[0] == "link":
            path.symlink_to(spec[1].replace("{arm}", str(arm)))
        elif spec[0] == "directory":
            path.mkdir()
        elif spec[0] == "mode":
            path.write_text(spec[2])
            path.chmod(spec[1])
    # A fixed staged modification time: only a copy that keeps it can match.
    for stage in arm.glob("stage*"):
        for path in stage.iterdir():
            os.utime(path, ns=(STAGED_NS, STAGED_NS))
    return arm


def script_run(producer, arm, arguments):
    """One setup step: the unchanged script, as the installer called it."""
    values = [value.replace("{arm}", str(arm)) for value in arguments]
    result = subprocess.run([sys.executable, producer, *values], cwd=arm, capture_output=True, text=True)
    if result.returncode != 0:
        raise SystemExit(f"setup step failed in {arm}: {values}: {result.stderr}")


def receipt_path(arm):
    return arm / "plugins/.priva-vep-plugins.json"


def edit_receipt(change):
    def step(arm, helper):
        value = json.loads(receipt_path(arm).read_text())
        change(value, arm)
        helper.write_receipt(receipt_path(arm), value)
    return step


def install_args(release, names=NAMES, transition=False, stage=None):
    values = ["install", "--directory", "{arm}/plugins", "--release", release, "--plugins", names,
              "--stage", "{arm}/" + (stage or f"stage{release}")]
    if transition:
        values += ["--transition-file", "{arm}/env/.priva-vep-release-transition.json"]
    return values


def utr_args(release="116"):
    return ["ensure-utr-reporting", "--directory", "{arm}/plugins", "--release", release]


UTR_NAMES = NAMES + "," + UTR
INSTALL_116 = ("run", install_args("116", UTR_NAMES))
ENSURE = ("run", utr_args())


def interrupted(stage):
    """Leave the transaction the script records before its live writes."""
    def step(arm, helper):
        after = json.loads(receipt_path(arm).read_text())
        before = json.loads(json.dumps(after))
        before.pop("local_patches")
        helper.write_receipt(arm / "plugins/.priva-utr-reporting.pending.json",
                             {"receipt_before": before, "receipt_after": after})
        if stage == "before_receipt":
            helper.write_receipt(receipt_path(arm), before)
    return step


def transition_config(**values):
    """The installer's configuration as JSON; "{arm}" is replaced per arm."""
    return json.dumps(values)


def catalogue():
    plugins = "{arm}/plugins"
    check = ["check", "--directory", plugins, "--release", "113", "--plugins", NAMES]
    transition_file = "{arm}/env/.priva-vep-release-transition.json"

    def transition(config, previous="113", release="116"):
        return ["transition", "--transition-file", transition_file, "--from-release", previous,
                "--release", release, "--config-json", config]

    installed_113 = {"files": stage_files("113"), "setup": [("run", install_args("113"))]}
    utr_ready = {"files": {**stage_files("116", utr=True), **stage_files("117", utr=True)},
                 "setup": [INSTALL_116]}
    utr_done = {**utr_ready, "setup": [INSTALL_116, ENSURE]}
    cases = {
        # Admission and main()'s own checks.
        "argument_unknown_mode": {"args": ["fetch", "--release", "116"]},
        "argument_missing_release": {"args": ["check", "--directory", plugins]},
        "release_not_numeric": {"args": ["check", "--directory", plugins, "--release", "v116"]},
        "release_leading_zero": {"args": ["check", "--directory", plugins, "--release", "0116"]},
        "directory_required": {"args": ["check", "--release", "116"]},
        "names_empty": {"args": ["check", "--directory", plugins, "--release", "116"]},
        "names_invalid": {"args": ["check", "--directory", plugins, "--release", "116", "--plugins", "A,1B"]},
        "directory_absent": {"args": ["check", "--directory", "{arm}/absent", "--release", "116",
                                      "--plugins", NAMES]},
        "stage_required": {"args": ["install", "--directory", plugins, "--release", "116", "--plugins", NAMES]},
        "transition_file_required": {"args": ["transition", "--release", "116"]},
        # check
        "check_without_receipt": {"args": check},
        "check_matching": {**installed_113, "args": check},
        "check_other_release": {**installed_113, "args": [*check[:4], "116", *check[5:]]},
        "check_changed_file": {**installed_113, "args": check,
                               "setup": [*installed_113["setup"], lambda arm, helper: (
                                   arm / "plugins/SpliceVault.pm").write_text(module("SpliceVault", " x"))]},
        "check_unrequested_plugin_missing": {**installed_113, "args": [*check[:6], NAMES + ",LOEUF"]},
        "check_schema_true": {**installed_113, "args": check,
                              "setup": [*installed_113["setup"],
                                        edit_receipt(lambda v, a: v.update(schema_version=True))]},
        "check_patches_null": {**installed_113, "args": check,
                               "setup": [*installed_113["setup"],
                                         edit_receipt(lambda v, a: v.update(local_patches=None))]},
        "check_patches_text": {**installed_113, "args": check,
                               "setup": [*installed_113["setup"],
                                         edit_receipt(lambda v, a: v.update(local_patches="SpliceAI.pm"))]},
        "check_config_digest_number": {**installed_113, "args": check,
                                       "setup": [*installed_113["setup"],
                                                 edit_receipt(lambda v, a: v.update(plugin_config_sha256=7))]},
        "check_recorded_utr_patch": {**utr_done, "args": ["check", "--directory", plugins, "--release", "116",
                                                          "--plugins", UTR_NAMES]},
        "check_utr_patch_backup_changed": {
            **utr_done, "args": ["check", "--directory", plugins, "--release", "116", "--plugins", UTR],
            "setup": [*utr_done["setup"], lambda arm, helper: (
                arm / "plugins/.priva-local-patches/UTRAnnotator116_all_effects/UTRAnnotator_all_effects.patch"
            ).write_text("changed\n")]},
        # install
        "install_fresh": {"files": {**stage_files("116"), "stage116/SpliceAI.pm": ("mode", 0o640, module("SpliceAI"))},
                          "args": install_args("116")},
        "install_duplicate_names": {"files": stage_files("116"), "args": install_args("116", "SpliceAI,SpliceAI")},
        "install_identical_unrecorded": {"files": {**stage_files("116"), "plugins/SpliceAI.pm": module("SpliceAI", " 116")},
                                         "args": install_args("116")},
        "install_customized_unrecorded": {"files": {**stage_files("116"), "plugins/SpliceAI.pm": module("SpliceAI", " x")},
                                          "args": install_args("116")},
        "install_destination_link": {"files": {**stage_files("116"), "plugins/SpliceAI.pm": ("link", "{arm}/absent")},
                                     "args": install_args("116")},
        "install_destination_directory": {"files": {**stage_files("116"), "plugins/SpliceAI.pm": ("directory",)},
                                          "args": install_args("116")},
        "install_wrong_package": {"files": {**stage_files("116"), "stage116/SpliceAI.pm": module("SpliceAI2")},
                                  "args": install_args("116")},
        "install_html_download": {"files": {**stage_files("116"), "stage116/SpliceVault.pm": "<html>busy</html>"},
                                  "args": install_args("116")},
        "install_package_spans_lines": {"files": {**stage_files("116"),
                                                  "stage116/SpliceAI.pm": "# x\n\n  package\n\tSpliceAI \n;\n1;\n"},
                                        "args": install_args("116")},
        "install_package_after_carriage_return": {"files": {**stage_files("116"),
                                                            "stage116/SpliceAI.pm": "# x\rpackage SpliceAI;\n"},
                                                  "args": install_args("116")},
        "install_not_utf8": {"files": {**stage_files("116"), "stage116/SpliceAI.pm": ("bytes", b"package SpliceAI;\xff\n")},
                             "args": install_args("116")},
        "install_stage_config_missing": {"files": {"stage116/SpliceAI.pm": module("SpliceAI"),
                                                   "stage116/SpliceVault.pm": module("SpliceVault")},
                                         "args": install_args("116")},
        "install_stage_plugin_empty": {"files": {**stage_files("116"), "stage116/SpliceVault.pm": ""},
                                       "args": install_args("116")},
        "install_release_switch": {**installed_113, "files": {**stage_files("113"), **stage_files("116")},
                                   "args": install_args("116")},
        "install_release_switch_customized": {
            **installed_113, "files": {**stage_files("113"), **stage_files("116")}, "args": install_args("116"),
            "setup": [*installed_113["setup"],
                      lambda arm, helper: (arm / "plugins/SpliceVault.pm").write_text(module("SpliceVault", " x"))]},
        "install_partial_retry": {
            **installed_113, "files": {**stage_files("113"), **stage_files("116")}, "args": install_args("116"),
            "setup": [*installed_113["setup"],
                      lambda arm, helper: shutil.copyfile(arm / "stage116/SpliceAI.pm", arm / "plugins/SpliceAI.pm")]},
        "install_switch_bad_metadata": {
            **installed_113, "files": {**stage_files("113"), **stage_files("116")}, "args": install_args("116"),
            "setup": [*installed_113["setup"], edit_receipt(lambda v, a: v.update(source="elsewhere"))]},
        "install_switch_old_file_missing": {
            **installed_113, "files": {**stage_files("113"), **stage_files("116")}, "args": install_args("116"),
            "setup": [*installed_113["setup"], lambda arm, helper: (arm / "plugins/SpliceVault.pm").unlink()]},
        "install_switch_unhashable_record": {
            **installed_113, "files": {**stage_files("113"), **stage_files("116")}, "args": install_args("116"),
            "setup": [*installed_113["setup"],
                      edit_receipt(lambda v, a: v["plugins"].update({"SpliceAI.pm": ["x"]}))]},
        "install_malformed_local_patches": {
            **installed_113, "files": {**stage_files("113")}, "args": install_args("113"),
            "setup": [*installed_113["setup"], edit_receipt(lambda v, a: v.update(local_patches=["x"]))]},
        "install_recorded_transition": {
            "files": {**stage_files("116"), "plugins/SpliceAI.pm": module("SpliceAI", " old")},
            "setup": [("run", transition(transition_config(vep_plugins_dir=plugins)))],
            "args": install_args("116", transition=True)},
        "install_transition_other_prefix": {
            "files": {**stage_files("116"), "plugins/SpliceAI.pm": module("SpliceAI", " old")},
            "setup": [("run", transition(transition_config(vep_plugins_dir=plugins))),
                      lambda arm, helper: os.replace(arm / "env/.priva-vep-release-transition.json",
                                                     arm / "plugins/moved.json")],
            "args": [*install_args("116"), "--transition-file", "{arm}/plugins/moved.json"]},
        "install_utr_keeps_patch": {**utr_done, "args": install_args("116", UTR_NAMES)},
        "install_utr_patch_needs_review": {**utr_done, "args": install_args("117", UTR_NAMES)},
        # transition
        "transition_plugins_dir": {
            "files": {"plugins/SpliceAI.pm": module("SpliceAI"), "plugins/lower_ok.pm": module("lower_ok"),
                      "plugins/.hidden.pm": "x", "plugins/1st.pm": "x", "plugins/notes.txt": "x",
                      "plugins/Linked.pm": ("link", "{arm}/plugins/SpliceAI.pm")},
            "args": transition(transition_config(vep_plugins_dir=plugins))},
        "transition_cache_root": {
            "files": {"cache/Custom/SpliceAI.pm": module("SpliceAI")},
            "args": transition(transition_config(cache_root="{arm}/x/../cache//",
                                                 cache_folders={"vep_plugin": "Custom"}))},
        "transition_cache_root_default_folder": {
            "files": {"cache/VEP_plugin/SpliceAI.pm": module("SpliceAI")},
            "args": transition(transition_config(cache_root="{arm}/cache", cache_folders=[]))},
        "transition_cache_root_relative": {"args": transition(json.dumps({"cache_root": "cache"}))},
        "transition_cache_root_slash": {"args": transition(json.dumps({"cache_root": "//"}))},
        "transition_cache_root_number": {"args": transition(json.dumps({"cache_root": 0}))},
        "transition_cache_root_empty": {"args": transition(transition_config(cache_root="",
                                                                             vep_plugins_dir=plugins))},
        "transition_folders_not_mapping": {"args": transition(json.dumps({"cache_root": "/c", "cache_folders": [1]}))},
        "transition_folder_parent": {"args": transition(json.dumps({"cache_root": "/c",
                                                                    "cache_folders": {"vep_plugin": ".."}}))},
        "transition_folder_null": {"args": transition(json.dumps({"cache_root": "/c",
                                                                  "cache_folders": {"vep_plugin": None}}))},
        "transition_same_release": {"args": transition("{}", previous="116")},
        "transition_previous_missing": {"args": ["transition", "--transition-file", transition_file,
                                                 "--release", "116", "--config-json", "{}"]},
        "transition_config_list": {"args": transition("[]")},
        "transition_config_invalid": {"args": transition("{")},
        "transition_config_empty_text": {"args": transition("")},
        "transition_plugins_dir_number": {"args": transition(json.dumps({"vep_plugins_dir": 5}))},
        "transition_plugins_dir_relative": {"args": transition(json.dumps({"vep_plugins_dir": "plugins"}))},
        "transition_directory_named_pm": {"files": {"plugins/Folder.pm": ("directory",)},
                                          "args": transition(transition_config(vep_plugins_dir=plugins))},
        "transition_prefix_missing": {"args": ["transition", "--transition-file", "{arm}/absent/t.json",
                                               "--from-release", "113", "--release", "116", "--config-json",
                                               json.dumps({"vep_plugins_dir": "/"})]},
        # clear-transition
        "clear_matching": {"files": {"plugins/SpliceAI.pm": module("SpliceAI")},
                           "setup": [("run", transition(transition_config(vep_plugins_dir=plugins)))],
                           "args": ["clear-transition", "--directory", plugins, "--release", "116",
                                    "--transition-file", transition_file]},
        "clear_other_release": {"setup": [("run", transition(transition_config(vep_plugins_dir=plugins)))],
                                "args": ["clear-transition", "--directory", plugins, "--release", "117",
                                         "--transition-file", transition_file]},
        "clear_without_file": {"args": ["clear-transition", "--directory", plugins, "--release", "116",
                                        "--transition-file", transition_file]},
        # ensure-utr-reporting
        "utr_not_applicable": {"args": utr_args("115")},
        "utr_directory_absent": {"args": ["ensure-utr-reporting", "--directory", "{arm}/absent", "--release", "116"]},
        "utr_install": {**utr_ready, "args": utr_args()},
        "utr_reused": {**utr_done, "args": utr_args()},
        "utr_without_receipt": {"files": {"plugins/UTRAnnotator.pm": ("copy", "tests/fixtures/UTRAnnotator116.original.pm")},
                                "args": utr_args()},
        "utr_target_changed": {**utr_ready, "args": utr_args(),
                               "setup": [INSTALL_116, lambda arm, helper: (arm / "plugins/UTRAnnotator.pm").write_text("x\n")]},
        "utr_unrecorded_customization": {
            **utr_ready, "args": utr_args(),
            "setup": [INSTALL_116, ENSURE, edit_receipt(lambda v, a: v.pop("local_patches"))]},
        "utr_pending_null": {**utr_done, "args": utr_args(),
                             "setup": [INSTALL_116, ENSURE, edit_receipt(lambda v, a: v.pop("local_patches")),
                                       lambda arm, helper: (arm / "plugins/.priva-utr-reporting.pending.json")
                                       .write_text("null\n")]},
        "utr_pending_invalid_json": {**utr_done, "args": utr_args(),
                                     "setup": [*utr_done["setup"],
                                               lambda arm, helper: (arm / "plugins/.priva-utr-reporting.pending.json")
                                               .write_text("{\n")]},
        "utr_pending_link": {**utr_done, "args": utr_args(),
                             "setup": [*utr_done["setup"], lambda arm, helper: (
                                 arm / "plugins/.priva-utr-reporting.pending.json").symlink_to(arm / "absent")]},
        "utr_resume_before_receipt": {**utr_done, "args": utr_args(),
                                      "setup": [*utr_done["setup"], interrupted("before_receipt")]},
        "utr_resume_after_receipt": {**utr_done, "args": utr_args(),
                                     "setup": [*utr_done["setup"], interrupted("after_receipt")]},
        "utr_resume_conflicting": {
            **utr_done, "args": utr_args(),
            "setup": [*utr_done["setup"], interrupted("after_receipt"),
                      edit_receipt(lambda v, a: v.update(note="changed"))]},
        "utr_backup_conflict": {
            **utr_ready, "args": utr_args(),
            "files": {**utr_ready["files"],
                      "plugins/.priva-local-patches/UTRAnnotator116_all_effects/UTRAnnotator_all_effects.patch": "x"}},
        "utr_backup_folder_link": {
            **utr_ready, "args": utr_args(),
            "files": {**utr_ready["files"], "elsewhere/x": "x",
                      "plugins/.priva-local-patches": ("link", "{arm}/elsewhere")}},
        "utr_recorded_other_profile": {
            **utr_done, "args": utr_args(),
            "setup": [*utr_done["setup"],
                      edit_receipt(lambda v, a: v["local_patches"]["UTRAnnotator.pm"].update(profile="other"))]},
        "utr_recorded_empty_patch": {
            **utr_ready, "args": utr_args(),
            "setup": [INSTALL_116, edit_receipt(lambda v, a: v.update(local_patches={"UTRAnnotator.pm": {}}))]},
        "utr_malformed_patches": {
            **utr_ready, "args": utr_args(),
            "setup": [INSTALL_116, edit_receipt(lambda v, a: v.update(local_patches={"UTRAnnotator.pm": 1}))]},
    }
    return cases


# ---------------------------------------------------------------- runs

def normalize(text, arm):
    if not isinstance(text, str):
        return text
    for spelling in sorted({str(arm), os.path.realpath(arm)}, key=len, reverse=True):
        text = text.replace(spelling, "<ARM>")
    return TEMPORARY.sub(r"\1<TEMP>", text)


def tree(root):
    state = {}
    for folder, directories, names in os.walk(root, followlinks=False):
        for name in sorted(directories + names):
            path = Path(folder) / name
            relative = normalize(str(path.relative_to(root)), root)
            if path.is_symlink():
                state[relative] = ["link", normalize(os.readlink(path), root)]
                continue
            mode = oct(path.stat().st_mode & 0o7777)
            if path.is_dir():
                state[relative] = ["directory", mode]
                continue
            data = path.read_bytes()
            try:
                content = ["text", normalize(data.decode("utf-8"), root)] if len(data) <= 65536 else None
            except UnicodeDecodeError:
                content = None
            content = content or ["bytes", hashlib.sha256(data).hexdigest()]
            entry = [*content, mode]
            if relative.startswith("plugins/") and relative.endswith(".pm") and "/" not in relative[8:]:
                staged = [stage / name for stage in sorted(root.glob("stage*")) if (stage / name).exists()]
                entry.append([stage.parent.name for stage in staged
                              if stage.stat().st_mtime_ns == path.stat().st_mtime_ns])
            state[relative] = entry
    return state


def run(command, folder, label, environment, cwd):
    with (folder / (label + ".stdout")).open("wb") as out, (folder / (label + ".stderr")).open("wb") as err:
        process = subprocess.run(list(map(str, command)), cwd=cwd, stdout=out, stderr=err,
                                 stdin=subprocess.DEVNULL, env=environment, timeout=600)
    return process.returncode


def outcome(folder, label, status, arm):
    stdout = (folder / (label + ".stdout")).read_bytes().decode("utf-8", "replace")
    stderr = (folder / (label + ".stderr")).read_bytes().decode("utf-8", "replace")
    result = {"exit": status, "stdout": normalize(stdout, arm), "tree": tree(arm)}
    if label == "python":
        recorded = folder / "exception.json"
        value = json.loads(recorded.read_text()) if recorded.exists() else {}
        if value:
            result.update(category=value.get("category"), message=value.get("message"))
    else:
        matched = re.search(r"(?:\A|\n)priva: ([A-Za-z_][A-Za-z_0-9]*): ((?:(?!\npriva: ).)*)\n\Z", stderr, re.S)
        if matched:
            result.update(category=matched[1], message=matched[2])
        elif stderr or status not in (0, 1, 2):
            result.update(category=None, message=stderr)
    if "message" in result:
        result["message"] = normalize(result["message"], arm)
    return result


def compare(python, native, before):
    if python["exit"] == 2 or native["exit"] == 2:
        # Argument admission: argparse usage text differs by design.
        return python["exit"] == native["exit"] == 2 and python["tree"] == native["tree"] == before
    left, right = dict(python), dict(native)
    left.pop("uncaught", None)
    return left == right


def check(args):
    producer = args.repo / "scripts/vep_plugins_release.py"
    if hashlib.sha256(producer.read_bytes()).hexdigest() != PRODUCER_SHA256:
        raise SystemExit("vep_plugins_release.py differs from the qualified source; update the pin")
    helper = load(producer, "qualified_vep_plugins_release")
    cases = catalogue()
    unknown = set(args.only or ()) - cases.keys()
    if unknown:
        raise SystemExit("Unknown cases: " + ", ".join(sorted(unknown)))
    args.work_root.mkdir(parents=True, exist_ok=False)
    rows = {}
    for name, case in cases.items():
        if args.only and name not in args.only:
            continue
        folder = args.work_root / name
        folder.mkdir()
        results, befores = {}, {}
        for label in ("python", "native"):
            arm = materialize((folder / label).absolute(), case, args.repo)
            for step in case.get("setup", []):
                if isinstance(step, tuple):
                    script_run(producer, arm, step[1])
                else:
                    step(arm, helper)
            befores[label] = tree(arm)
            environment = {key: value for key, value in os.environ.items() if not key.startswith("PRIVA_")}
            environment["PYTHONDONTWRITEBYTECODE"] = "1"
            arguments = [value.replace("{arm}", str(arm)) for value in case["args"]]
            if label == "python":
                command = [sys.executable, Path(__file__).resolve(), "oracle", "--producer", producer,
                           "--exception", folder / "exception.json", "--", *arguments]
            else:
                environment["PRIVA_DATA_DIR"] = str(args.repo / "data")
                command = [args.binary, "vep-plugins", *arguments]
            status = run(command, folder, label, environment, arm)
            results[label] = outcome(folder, label, status, arm)
        passed = befores["python"] == befores["native"] and \
            compare(results["python"], results["native"], befores["native"])
        rows[name] = {"passed": passed, "python": results["python"], "native": results["native"]}
        write_json(folder / "comparison.json", rows[name])
    failed = sorted(name for name, row in rows.items() if not row["passed"])
    report = {
        "producer": str(producer), "producer_sha256": PRODUCER_SHA256,
        "binary": str(args.binary),
        "binary_sha256": hashlib.sha256(Path(args.binary).read_bytes()).hexdigest(),
        "cases": len(rows), "python_successes": sum(row["python"]["exit"] == 0 for row in rows.values()),
        "failed": failed, "rows": rows,
    }
    write_json(args.work_root / "complete.json", report)
    print(json.dumps({key: report[key] for key in ("cases", "python_successes", "failed")}))
    if failed:
        raise SystemExit(1)


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    commands = parser.add_subparsers(dest="procedure", required=True)
    p = commands.add_parser("check")
    p.add_argument("--binary", type=Path, required=True)
    p.add_argument("--repo", type=Path, required=True)
    p.add_argument("--work-root", type=Path, required=True)
    p.add_argument("--only", action="append")
    p.set_defaults(function=check)
    p = commands.add_parser("oracle")
    p.add_argument("--producer", type=Path, required=True)
    p.add_argument("--exception", type=Path, required=True)
    p.add_argument("arguments", nargs=argparse.REMAINDER)
    p.set_defaults(function=oracle)
    args = parser.parse_args()
    if args.procedure == "oracle" and args.arguments[:1] == ["--"]:
        args.arguments = args.arguments[1:]
    args.function(args)


if __name__ == "__main__":
    main()
