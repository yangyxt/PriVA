#!/usr/bin/env python3
"""Runs one install step both ways and compares what each records.

Old: `priva install-task run` with the registry before the step became its
own process (tests/fixtures/install-steps-1879c0c.json), which calls the
step's shell function from scripts/install_utils.sh. New: the step's process
(conf/install_steps.json "module"), run by Nextflow on the same plan. Both get
the same configuration, the same receipts of the steps it needs (--input) and,
with --net, the same server: tests/fake_net.py answers curl and wget from that
folder (HOST/PATH), so a step can be compared without the network.

The output keys must be the same, and each key's file the same bytes (a
folder: the same files with the same bytes, git's own folder and the commit
file the process records left out; a gzip file that is not BGZF may differ in
its header, which holds the name and time it was written with, and then its
decompressed bytes must be the same). Exit 1 on any difference.

  python3 tests/install_step_compare.py --step B6 --config C.yaml --work W \\
      --priva rust/target/release/priva --nextflow nextflow [--net DIR] \\
      [--input RECEIPT]... [--new-only]

--new-only runs only the step's process and prints the keys it records (for a
case the old function handled by writing beside a file you supply).

A key whose bytes differ can be compared again apart from what is bound to
differ between two runs, and is then reported as the same apart from it:

--volatile NAME (repeatable): without the members named NAME at any depth of
a JSON file (or gzip JSON), such as the time a builder stamps on what it
writes (built_at), and without the column NAME of a tab-separated file whose
first line names it (a download log's downloaded_at).
--own-paths: with each route's own step folder (the receipt's "folder") read
as the same, for a file that names where it lives.
A VCF's meta lines are volatile by key: --volatile 'bcftools_*Command' leaves
out the ##bcftools_viewCommand=... lines (the paths and time of the command
that wrote it), compared as shell patterns. Such a file's index then differs
with it, and is only checked to be there on both sides.

--family: a key that names the chrX member of a chromosome family
(NAME.chrX.vcf.gz and the like) is compared with its 23 other members.
"""
import argparse
import fnmatch
import gzip
import hashlib
import json
import os
import shutil
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
OLD_REGISTRY = ROOT / "tests" / "fixtures" / "install-steps-1879c0c.json"
LINKED = ("conf", "scripts", "data", "rust", "modules", "workflows", "bin", "cache_resources.json")


def sha(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for block in iter(lambda: f.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


def fingerprint(path: Path) -> str:
    """A file's SHA-256, or a folder's files and their SHA-256 (no .git)."""
    if path.is_file():
        return sha(path)
    lines = []
    for item in sorted(path.rglob("*")):
        rel = item.relative_to(path)
        if rel.parts and rel.parts[0] == ".git":
            continue
        if item.is_file():
            lines.append(f"{rel} {sha(item)}")
    return hashlib.sha256("\n".join(lines).encode()).hexdigest() + f" ({len(lines)} files)"


def gunzipped(path: Path) -> str | None:
    """The SHA-256 of a gzip file's decompressed bytes (None if not gzip)."""
    if not path.is_file() or path.read_bytes()[:2] != b"\x1f\x8b":
        return None
    h = hashlib.sha256()
    with gzip.open(path, "rb") as f:
        for block in iter(lambda: f.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


def without(value, volatile):
    """A JSON value without the members named in volatile, at any depth."""
    if isinstance(value, dict):
        return {k: without(v, volatile) for k, v in value.items() if k not in volatile}
    if isinstance(value, list):
        return [without(v, volatile) for v in value]
    return value


def normalized(path: Path, folder: str, volatile) -> bytes:
    """A file's bytes (decompressed) with its folder read as <FOLDER>, and
    without its volatile JSON members or tab-separated columns."""
    raw = path.read_bytes()
    if raw[:2] == b"\x1f\x8b":
        raw = gzip.decompress(raw)
    if folder:
        raw = raw.replace(folder.encode(), b"<FOLDER>")
    if not volatile:
        return raw
    if raw.startswith(b"##fileformat=VCF"):
        def kept(line: bytes) -> bool:
            if not line.startswith(b"##"):
                return True
            key = line[2:].split(b"=", 1)[0].decode(errors="replace")
            return not any(fnmatch.fnmatchcase(key, pattern) for pattern in volatile)
        return b"\n".join(line for line in raw.split(b"\n") if kept(line))
    try:
        value = json.loads(raw)
    except ValueError:
        value = None
    if isinstance(value, (dict, list)):
        return json.dumps(without(value, set(volatile)), separators=(",", ":")).encode()
    lines = raw.split(b"\n")
    header = lines[0].split(b"\t")
    drop = {i for i, name in enumerate(header) if name.decode(errors="replace") in volatile}
    if not drop:
        return raw
    return b"\n".join(b"\t".join(c for i, c in enumerate(line.split(b"\t")) if i not in drop) for line in lines)


def apart(path: Path, folder: str, volatile) -> str:
    """The fingerprint of a file or folder, compared apart (see --volatile)."""
    files = [path] if path.is_file() else sorted(
        item for item in path.rglob("*")
        if item.is_file() and not (item.relative_to(path).parts and item.relative_to(path).parts[0] == ".git"))
    lines = [f"{item.relative_to(path) if item != path else ''} "
             f"{hashlib.sha256(normalized(item, folder, volatile)).hexdigest()}" for item in files]
    return hashlib.sha256("\n".join(lines).encode()).hexdigest()


def run(command, **kw):
    print("+", " ".join(str(c) for c in command), file=sys.stderr)
    return subprocess.run([str(c) for c in command], **kw)


def plan(args, work: Path) -> Path:
    releases = work / "releases"
    releases.mkdir(parents=True, exist_ok=True)
    registry = json.loads((ROOT / "conf" / "install_steps.json").read_text())
    files = []
    for step in registry["steps"]:
        if step.get("release"):
            out = releases / f"release_{step['id']}.json"
            run([args.priva, "install-task", "release", "--source-root", ROOT, "--registry",
                 ROOT / "conf" / "install_steps.json", "--config", args.config, "--step", step["id"],
                 "--out", out, "--stub"], check=True)
            files += ["--releases", out]
    run([args.priva, "install-task", "plan", "--source-root", ROOT, "--registry",
         ROOT / "conf" / "install_steps.json", "--config", args.config, "--out", work / "plan",
         "--refresh", ""] + files, check=True)
    return work / "plan" / "plan.json"


def environment(args) -> dict:
    env = dict(os.environ)
    if args.net:
        net = Path(args.net).resolve()
        subprocess.run([sys.executable, Path(__file__).with_name("fake_net.py"), "install", net], check=True)
        env["PATH"] = f"{net / 'bin'}:{env['PATH']}"
        env["FAKE_NET_ROOT"] = str(net)
    return env


def old_route(args, plan_json: Path, work: Path, env: dict) -> dict:
    here = work / "old" / "run"
    here.mkdir(parents=True, exist_ok=True)
    inputs = [x for r in args.input for x in ("--input", Path(r).resolve())]
    run([args.priva, "install-task", "run", "--source-root", ROOT, "--registry", OLD_REGISTRY,
         "--plan", plan_json, "--step", args.step, "--lineage", "old", "--folder", work / "old" / "folder",
         "--receipt", "receipt.json"] + inputs,
        cwd=here, env=dict(env, PRIVA_BINARY=str(Path(args.priva).resolve()), INSTALL_THREADS="2"), check=True)
    return json.loads((here / "receipt.json").read_text())


def new_route(args, plan_json: Path, work: Path, env: dict) -> dict:
    here = work / "new"
    here.mkdir(parents=True, exist_ok=True)
    for name in LINKED:
        link = here / name
        if not link.exists():
            link.symlink_to(ROOT / name)
    registry = json.loads((ROOT / "conf" / "install_steps.json").read_text())
    entry = next(s for s in registry["steps"] if s["id"] == args.step)
    module = ROOT / entry["module"]
    stage = ", ".join(f"'{k}'" for k in entry.get("stage", []))
    (here / "harness.nf").write_text(f"""
include {{ installNative }} from '{ROOT}/workflows/install'
include {{ {entry.get('process', entry['name'])} as {entry['name']} }} from '{str(module)[:-3]}'
workflow {{
    def plan = Channel.of(file(params.plan))
    def needs = (params.inputs ?: '').tokenize(',').collect {{ r -> Channel.of(file(r)) }}
    {entry['name']}(installNative('{args.step}', plan, needs, [{stage}]))
}}
""")
    (here / "limits.config").write_text("process.resourceLimits = [cpus: 4, memory: 12.GB, time: 2.h]\n")
    resource_root = here / "res"
    inputs = ",".join(str(Path(r).resolve()) for r in args.input)
    command = [args.nextflow, "run", "harness.nf", "-c", "limits.config", "-ansi-log", "false",
               "--plan", plan_json, "--resource_root", resource_root,
               "--priva_binary", Path(args.priva).resolve()]
    if inputs:
        command += ["--inputs", inputs]
    result = run(command, cwd=here, env=dict(env, NXF_SYNTAX_PARSER="v2", NXF_OFFLINE="true"))
    if result.returncode != 0:
        raise SystemExit(f"the {args.step} process failed; see {here}/.nextflow.log and work/")
    receipts = sorted((resource_root / args.step).glob("*/receipt.json"))
    return json.loads(receipts[-1].read_text())


def compare(args, old: dict, new: dict, key: str, a: str, b: str) -> int:
    """Compare one key's files of the two routes; 1 when they differ."""
    if not (Path(a).exists() and Path(b).exists()):
        print(f"  {key}: missing on {'both sides' if not Path(a).exists() and not Path(b).exists() else 'one side'}")
        print(f"    old {a}\n    new {b}")
        return 1
    fa, fb = fingerprint(Path(a)), fingerprint(Path(b))
    same = fa == fb
    if not same and gunzipped(Path(a)) is not None and gunzipped(Path(a)) == gunzipped(Path(b)):
        # A gzip header holds the name and time it was written with, which
        # differ between two runs of the same command.
        print(f"  {key}: same decompressed {gunzipped(Path(a))} (gzip header differs)")
        print(f"    old {a}\n    new {b}")
        return 0
    apart_from = False
    if not same and (args.volatile or args.own_paths):
        folders = (old.get("folder", ""), new.get("folder", "")) if args.own_paths else ("", "")
        pa, pb = apart(Path(a), folders[0], args.volatile), apart(Path(b), folders[1], args.volatile)
        if pa == pb:
            reasons = args.volatile + (["its own folder"] if args.own_paths else [])
            print(f"  {key}: same apart from {', '.join(reasons)} {pa}")
            print(f"    old {a}\n    new {b}")
            same = apart_from = True
    if not apart_from:
        print(f"  {key}: {'same' if same else 'DIFFERENT'} {fa}" + ("" if same else f" / {fb}"))
        print(f"    old {a}\n    new {b}")
    differ = not same
    # The index files beside a key's file are read with it.
    for ext in (".tbi", ".csi", ".fai", ".gzi"):
        ia, ib = Path(a + ext), Path(b + ext)
        if ia.exists() or ib.exists():
            ok = ia.exists() and ib.exists() and (apart_from or sha(ia) == sha(ib))
            differ += not ok
            what = "there on both sides" if apart_from else "same"
            print(f"  {key}{ext}: {what if ok else 'DIFFERENT or missing on one side'}")
    return int(bool(differ))


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--step", required=True)
    parser.add_argument("--config", required=True)
    parser.add_argument("--work", required=True)
    parser.add_argument("--priva", required=True)
    parser.add_argument("--nextflow", required=True)
    parser.add_argument("--net")
    parser.add_argument("--input", action="append", default=[])
    parser.add_argument("--new-only", action="store_true")
    parser.add_argument("--volatile", action="append", default=[])
    parser.add_argument("--own-paths", action="store_true")
    parser.add_argument("--family", action="store_true")
    args = parser.parse_args()
    args.config = str(Path(args.config).resolve())
    args.priva = str(Path(args.priva).resolve())
    work = Path(args.work).resolve()
    if work.exists():
        subprocess.run(["chmod", "-R", "u+w", work])
        shutil.rmtree(work)
    work.mkdir(parents=True)
    env = environment(args)
    plan_json = plan(args, work)
    if args.new_only:
        new = new_route(args, plan_json, work, env)
        for key, value in sorted(new["keys"].items()):
            print(f"  {key}: {fingerprint(Path(value))} {value}")
        return
    old = old_route(args, plan_json, work, env)
    new = new_route(args, plan_json, work, env)
    differ = 0
    keys = sorted(set(old["keys"]) | set(new["keys"]))
    print(f"{args.step}: {len(keys)} key(s)")
    for key in keys:
        a, b = old["keys"].get(key), new["keys"].get(key)
        if a is None or b is None:
            print(f"  {key}: only {'new' if a is None else 'old'} records it")
            differ += 1
            continue
        pairs = [(key, a, b)]
        if args.family and ".chrX." in Path(a).name and ".chrX." in Path(b).name:
            for chrom in [str(n) for n in range(1, 23)] + ["Y"]:
                member = lambda path: str(Path(path).with_name(Path(path).name.replace(".chrX.", f".chr{chrom}.")))
                pairs.append((f"{key}[chr{chrom}]", member(a), member(b)))
        for label, a, b in pairs:
            differ += compare(args, old, new, label, a, b)
    sys.exit(1 if differ else 0)


if __name__ == "__main__":
    main()
