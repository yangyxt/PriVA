"""A stand-in priva program for installer tests that have no Rust build.

The installer runs native builders through PRIVA_BINARY. This program prints a
priva version, appends each call's arguments as one JSON line to
PRIVA_FAKE_CALLS, and then runs the Python script the native command replaced,
with the equivalent arguments, through `python` in a new Bash. A test's own
exported `python` function (export -f python) therefore still sees the call.
PRIVA_FAKE_STATUS makes every builder call exit with that status instead.
`builder-identity NAME` prints a fixed identity for NAME, or
PRIVA_FAKE_BUILDER_SHA256 when it is set, and runs no script.
This program is never used by an installer.
"""
import json
from pathlib import Path
import re
import sys
import tempfile

ROOT = Path(__file__).resolve().parents[1]

PROGRAM = r'''#!PYTHON
import hashlib, json, os, subprocess, sys
SCRIPTS = SCRIPT_FOLDER
args = sys.argv[1:]
if args == ["--version"]:
    print("priva fixture")
    raise SystemExit(0)
if os.environ.get("PRIVA_FAKE_CALLS"):
    with open(os.environ["PRIVA_FAKE_CALLS"], "a") as calls:
        calls.write(json.dumps(args) + "\n")
if os.environ.get("PRIVA_FAKE_STATUS"):
    raise SystemExit(int(os.environ["PRIVA_FAKE_STATUS"]))


def options(values, names):
    """Return the --name value pairs, refusing any name the command lacks."""
    if len(values) % 2:
        raise SystemExit("priva fixture: odd option list " + repr(values))
    pairs = list(zip(values[::2], values[1::2]))
    for name, _ in pairs:
        if name not in names:
            raise SystemExit("priva fixture: unknown option " + name)
    return pairs


command, rest = args[0], args[1:]
environment = dict(os.environ)
if command == "builder-identity" and len(rest) == 1:
    # The native value covers the builder's sources; the fixture names one per builder.
    value = os.environ.get("PRIVA_FAKE_BUILDER_SHA256") or hashlib.sha256(rest[0].encode()).hexdigest()
    print(json.dumps({"builder": rest[0], "builder_source_sha256": value}, separators=(",", ":")))
    raise SystemExit(0)
if command == "installed-environment":
    script, forwarded = "installed_environment.py", rest
elif command == "indexed-plugin-cache":
    flags = [value for value in rest[1:] if value == "--print-path"]
    pairs = options([value for value in rest[1:] if value != "--print-path"],
                    {"--input", "--kind", "--cache-dir", "--assembly", "--threads", "--tmp-dir"})
    script, forwarded = "indexed_plugin_cache.py", [rest[0], *sum(pairs, ())] + flags
elif command == "splicevault-convert":
    pairs = options(rest, {"--input_tsv", "--output_tsv", "--chain_file", "--reference_fasta",
                           "--temp_dir", "--chunk_size"})
    script, forwarded = "SpliceVault_tsv_vcf_conversion.py", list(sum(pairs, ()))
elif command == "conservation-resource" and rest[:1] in (["validate"], ["publish"]):
    pairs = options(rest[1:], {"--input", "--assembly", "--reference", "--tool-prefix",
                               "--source-url", "--publish-to"})
    if "--tool-prefix" not in dict(pairs) or ("--publish-to" in dict(pairs)) != (rest[0] == "publish"):
        raise SystemExit("priva fixture: conservation-resource procedure and options disagree")
    script = "conservation_resource_check.py"
    forwarded = list(sum((pair for pair in pairs if pair[0] != "--tool-prefix"), ()))
elif command == "hpo-assertions" and rest[:1] == ["build"]:
    values = dict(options(rest[1:], {"--genes-to-phenotype", "--phenotype-hpoa", "--output",
                                     "--minimum-match-rate", "--disease-scope-registry"}))
    script = "collapse_HPO_anno.py"
    forwarded = [values.pop("--genes-to-phenotype"), values.pop("--phenotype-hpoa"),
                 values.pop("--output"), *sum(values.items(), ())]
elif command == "mondo-scope" and rest[:1] == ["build"]:
    pairs = options(rest[1:], {"--mondo-obo", "--hpo-assertions", "--registry-output",
                               "--annotated-assertions-output", "--manual-overrides"})
    script, forwarded = "build_mondo_disease_scope.py", list(sum(pairs, ()))
elif command == "clingen-map" and rest[:1] == ["build"]:
    short = {"--clinvar-vcf": "-v", "--clingen-tab": "-t", "--output": "-o"}
    pairs = options(rest[1:], set(short))
    script = "prepare_clingen_map.py"
    forwarded = list(sum(((short[name], value) for name, value in pairs), ()))
elif command == "clinvar-contract" and rest[:1] == ["header"]:
    values = dict(options(rest[1:], {"--vcf", "--stage"}))
    script = "clinvar_annotation_contract.py"
    forwarded = [values["--vcf"], "--stage", values.get("--stage", "complete")]
elif command == "clinvar-contract" and rest[:1] == ["cache"]:
    values = dict(options(rest[1:], {"--input", "--kind"}))
    script = "clinvar_annotation_contract.py"
    forwarded = [values["--input"], "--cache-kind", values["--kind"]]
elif command == "clinvar-contract" and rest[:1] == ["af-paths-check"]:
    pairs = options(rest[1:], {"--vcf", "--config", "--gene", "--exon"})
    script, forwarded = "clinvar_annotation_contract.py", ["af-paths-check", *sum(pairs, ())]
elif command == "clinvar-gene-af" and rest[:1] == ["build"]:
    values = dict(options(rest[1:], {"--vcf", "--output", "--threads"}))
    script = "stat_gene_patho_afs.py"
    forwarded = [values["--vcf"], values["--output"], values["--threads"]]
elif command == "clinvar-exon-af" and rest[:1] == ["build"]:
    long = {"--vcf": "--vcf_path", "--output": "--output", "--threads": "--threads"}
    pairs = options(rest[1:], set(long))
    script = "stat_exon_patho_afs.py"
    forwarded = list(sum(((long[name], value) for name, value in pairs), ()))
elif command == "clinvar-aa-splice" and rest[:1] == ["build"]:
    values = dict(options(rest[1:], {"--vcf", "--aa-output", "--splice-output", "--threads"}))
    script = "stat_aachange_clinvar.py"
    forwarded = [values["--vcf"], values["--aa-output"], values["--splice-output"], values["--threads"]]
elif command == "clinvar-gene-stats" and rest[:1] == ["build"]:
    values = dict(options(rest[1:], {"--vcf", "--output", "--python-hash-seed", "--clinvar-summary",
                                     "--mane-summary", "--threads", "--config"}))
    # The native seed names the Python run whose set order it reproduces.
    environment["PYTHONHASHSEED"] = values.pop("--python-hash-seed")
    long = {"--vcf": "--clinvar_anno_vcf", "--output": "--output_pickle"}
    script = "clinvar_stat_variants.py"
    forwarded = list(sum(((long.get(name, name), value) for name, value in values.items()), ()))
else:
    raise SystemExit("priva fixture: no Python equivalent for " + repr(args))
result = subprocess.run(["bash", "--noprofile", "--norc", "-c", 'python "$@"', "priva-fixture",
                         os.path.join(SCRIPTS, script), *forwarded], env=environment)
raise SystemExit(result.returncode)
'''


def install(directory=None, scripts=ROOT / "scripts"):
    """Write the program into DIRECTORY (a new temporary folder by default)."""
    folder = Path(directory or tempfile.mkdtemp(prefix="priva_fake_binary_"))
    program = folder / "priva"
    program.write_text(PROGRAM.replace("#!PYTHON", "#!" + sys.executable, 1)
                       .replace("SCRIPT_FOLDER", repr(str(scripts)), 1))
    program.chmod(0o755)
    return program


def calls(path):
    """Return the recorded argument lists, oldest first."""
    path = Path(path)
    if not path.exists():
        return []
    return [json.loads(line) for line in path.read_text().splitlines()]


def environment(directory, scripts=ROOT / "scripts", **extra):
    """Return PRIVA_BINARY and PRIVA_FAKE_CALLS settings for one test folder.

    priva_resolve_binary accepts only paths without spaces or shell characters,
    so a folder name outside that set gets its program in a new temporary folder.
    """
    folder = Path(directory)
    home = folder / "fake_priva"
    if not re.fullmatch(r"[A-Za-z0-9_./:+-]+", str(home.absolute())):
        home = Path(tempfile.mkdtemp(prefix="priva_fake_binary_"))
    home.mkdir(exist_ok=True)
    return dict(PRIVA_BINARY=str(install(home, scripts).absolute()),
                PRIVA_FAKE_CALLS=str(folder / "priva_calls.jsonl"), **extra)
