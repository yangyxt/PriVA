"""A stand-in priva program for installer tests that have no Rust build.

The installer runs native builders through PRIVA_BINARY. This program prints a
priva version, appends each call's arguments as one JSON line to
PRIVA_FAKE_CALLS, and then runs the Python script the native command replaced,
with the equivalent arguments, through `python` in a new Bash. A test's own
exported `python` function (export -f python) therefore still sees the call.
PRIVA_FAKE_STATUS makes every builder call exit with that status instead.
`builder-identity NAME` prints a fixed identity for NAME, or
PRIVA_FAKE_BUILDER_SHA256 when it is set, and runs no script.
`cache-install resolve` prints the build config's value for the key (relative
values under --source-root) without converting it, and runs no script.
`mapped-artifact run/recover --native-producer pm1-regions` runs the Python
publication helper with the Python PM1 producer and its --skip_das arguments.
`das-install-state` runs scripts/das_install_state.py with the same arguments.
`installer-prerequisites` runs scripts/installer_prerequisites.py with the same
arguments (CONFIG first, then the script's own options).
`reference-resources` runs scripts/reference_resources.py with the same
subcommand, CONFIG and options.
`vep-plugins MODE` runs scripts/vep_plugins_release.py with the same mode and
options; `gnomad-cache COMMAND` runs scripts/gnomad_cache.py and
`alphamissense-cache COMMAND` runs scripts/alphamissense_cache.py and
`spliceai2-prepare COMMAND` runs scripts/spliceai2_to_vep_vcf.py and
`clinvar-gene-refresh COMMAND` runs scripts/clinvar_gene_refresh.py likewise.
`gene-nonlof-mechanism` runs scripts/build_gene_nonlof_mechanism_cache.py and
`hpo-condition-mechanism` runs scripts/build_hpo_condition_mechanism_cache.py
with every option except --source-root.
`installer-config-paths resolve/managed-roots` runs scripts/installer_config_paths.py
with the same subcommand and options; it needs an explicit PYTHONHASHSEED, as the
native command does.
The AlphaMissense statistics, exon-history, domain-metadata and das-domains
transcripts/analyze commands replace no single Python script: the fixture
writes small JSON stand-ins that record the inputs and the current builder
identity, and their check/identity/export/merge commands read them back.
PRIVA_FAKE_FAIL="COMMAND SUBCOMMAND" makes that one call exit 73 before writing;
PRIVA_FAKE_CHANGE=PATH makes `das-domains analyze` append to PATH first.
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
if os.environ.get("PRIVA_FAKE_FAIL") == " ".join(args[:2]):
    raise SystemExit(73)


def current(name):
    """The fixture builder identity, as builder-identity prints it."""
    return os.environ.get("PRIVA_FAKE_BUILDER_SHA256") or hashlib.sha256(name.encode()).hexdigest()


def stand_in(path, **value):
    with open(path, "w") as stream:
        json.dump(value, stream, sort_keys=True)


def read(path):
    with open(path) as stream:
        return json.load(stream)


SYNTHETIC = {"alphamissense-domain": "alphamissense-domain", "domain-metadata": "alphamissense-domain",
             "alpha-exon-history": "alpha-exon-history"}
if command in SYNTHETIC and rest[:1] in (["chromosome"], ["merge"], ["build"], ["check"], ["identity"],
                                         ["export"]):
    builder = current(SYNTHETIC[command])
    procedure, pairs = rest[0], options(rest[1:], {"--vcf", "--chrom", "--output", "--input", "--hash-seed",
                                                    "--builder-sha256", "--assembly", "--threads"})
    values = dict(pairs)
    inputs = [value for name, value in pairs if name == "--input"]
    if procedure == "chromosome":
        stand_in(values["--output"], command=command, stage="raw", vcf=values["--vcf"],
                 chrom=values["--chrom"], builder_source_sha256=builder)
    elif procedure == "build":
        receipt = values["--output"] + ".source.json"
        if os.path.exists(values["--output"]) and os.path.exists(receipt) and read(receipt) == values:
            print(json.dumps({"reused": True}))
            raise SystemExit(0)
        stand_in(values["--output"], command=command, stage="final", vcf=values["--vcf"],
                 builder_source_sha256=builder)
        stand_in(receipt, **values)
    elif procedure == "merge":
        shards = [read(path) for path in inputs]
        if any(shard["stage"] != "raw" or shard["builder_source_sha256"] != builder for shard in shards):
            raise SystemExit("priva fixture: checkpoint producer/implementation differs")
        stand_in(values["--output"], command=command, stage="final", vcf=shards[0]["vcf"],
                 chromosomes=[shard["chrom"] for shard in shards], seed=values.get("--hash-seed"),
                 builder_source_sha256=builder)
    elif procedure == "export":
        recorded = read(inputs[0])
        if recorded["builder_source_sha256"] != builder or recorded["stage"] != "final":
            raise SystemExit("priva fixture: only final current history can be exported")
        with open(values["--output"], "w") as stream:
            stream.write("# " + json.dumps(recorded, sort_keys=True) + "\n")
    else:
        recorded = read(inputs[0])
        if procedure == "identity":
            print(json.dumps({"builder_source_sha256": recorded["builder_source_sha256"],
                              "kind": command, "stage": recorded["stage"]}, separators=(",", ":")))
            raise SystemExit(0)
        if "--vcf" in values and values["--vcf"] != recorded["vcf"]:
            raise SystemExit("priva fixture: selected VCF/index association differs")
        if recorded["builder_source_sha256"] != values.get("--builder-sha256", builder):
            raise SystemExit("priva fixture: selected producer association differs")
        print(json.dumps({"passed": True}))
    raise SystemExit(0)
if command == "das-domains" and rest[:1] == ["transcripts"]:
    values = dict(options(rest[1:], {"--statistics", "--exon-history", "--source-builder-sha256",
                                     "--python-hash-seed", "--output", "--history-builder-sha256"}))
    stand_in(values["--output"], builder_source_sha256=current("das-domains"), **values)
    raise SystemExit(0)
if command == "das-domains" and rest[:1] == ["analyze"]:
    values = dict(options(rest[1:], {"--statistics", "--source-builder-sha256", "--interpro",
                                     "--normalization-json", "--python-hash-seed", "--threads", "--assembly",
                                     "--output-dir", "--metadata", "--spots", "--clinvar-vcf"}))
    if os.environ.get("PRIVA_FAKE_CHANGE"):
        with open(os.environ["PRIVA_FAKE_CHANGE"], "ab") as stream:
            stream.write(b"changed during analysis")
    for stem, suffix in (("domain_tolerance_analysis", "tsv"), ("domain_filter_results", "tsv"),
                         ("curated_functional_references", "json"), ("excluded_domains", "json"),
                         ("all_intolerant_domains", "rkyv")):
        stand_in(os.path.join(values["--output-dir"], stem + "." + values["--assembly"] + "." + suffix),
                 output=stem, builder_source_sha256=current("das-domains"), **values)
    raise SystemExit(0)
if command == "cache-install" and rest[:1] == ["resolve"]:
    values = dict(options(rest[1:], {"--key", "--build-config", "--source-root"}))
    text = open(values["--build-config"]).read()
    try:
        document = json.loads(text)
    except ValueError:
        import yaml
        document = yaml.safe_load(text)
    selected = document.get(values["--key"])
    if not selected:
        raise SystemExit("priva fixture: " + values["--key"] + " is not configured")
    print(os.path.join(values["--source-root"], selected))
    raise SystemExit(0)
if command == "installed-environment":
    script, forwarded = "installed_environment.py", rest
elif command == "das-install-state":
    script, forwarded = "das_install_state.py", rest
elif command == "installer-prerequisites":
    # Same positional CONFIG and options; the native command requires --source-root.
    if "--source-root" not in rest:
        raise SystemExit("priva fixture: installer-prerequisites needs --source-root")
    script, forwarded = "installer_prerequisites.py", rest
elif command == "reference-resources" and rest[:1] in (["sources"], ["preflight"], ["install"]):
    # Same subcommand, CONFIG and options; the native command requires --data-dir
    # for preflight and install.
    if rest[0] != "sources" and "--data-dir" not in rest:
        raise SystemExit("priva fixture: reference-resources needs --data-dir")
    script, forwarded = "reference_resources.py", rest
elif command == "installer-config-paths" and rest[:1] in (["resolve"], ["managed-roots"]):
    flags = [value for value in rest[1:] if value == "--check"]
    names = {"--config", "--source-root"}
    if rest[0] == "resolve":
        names |= {"--managed-updates-json", "--output"}
    pairs = options([value for value in rest[1:] if value != "--check"], names)
    if "--source-root" not in dict(pairs) or not os.environ.get("PYTHONHASHSEED", "").isdigit():
        raise SystemExit("priva fixture: installer-config-paths needs --source-root and PYTHONHASHSEED")
    script, forwarded = "installer_config_paths.py", [rest[0], *sum(pairs, ()), *flags]
elif command == "vep-plugins" and rest[:1] in (["check"], ["install"], ["transition"],
                                              ["clear-transition"], ["ensure-utr-reporting"]):
    script, forwarded = "vep_plugins_release.py", rest
elif command == "gnomad-cache" and rest[:1] in (["members"], ["validate"], ["default-path"]):
    script, forwarded = "gnomad_cache.py", rest
elif command == "alphamissense-cache" and rest[:1] in (["install"], ["validate"], ["check-annotation"]):
    script, forwarded = "alphamissense_cache.py", rest
elif command == "spliceai2-prepare" and rest[:1] in (["prepare"], ["gene-context"]):
    script, forwarded = "spliceai2_to_vep_vcf.py", rest
elif command == "clinvar-gene-refresh":
    script, forwarded = "clinvar_gene_refresh.py", rest
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
elif command in ("gene-nonlof-mechanism", "hpo-condition-mechanism"):
    # Every option but --source-root passes through unchanged, flags included.
    script = {"gene-nonlof-mechanism": "build_gene_nonlof_mechanism_cache.py",
              "hpo-condition-mechanism": "build_hpo_condition_mechanism_cache.py"}[command]
    forwarded = []
    remaining = list(rest)
    while remaining:
        token = remaining.pop(0)
        if token == "--source-root":
            remaining.pop(0)
        else:
            forwarded.append(token)
elif command == "clinvar-contract" and rest[:1] == ["legacy-cache"]:
    values = dict(options(rest[1:], {"--input", "--kind"}))
    script = "clinvar_annotation_contract.py"
    forwarded = [values["--input"], "--cache-kind", values["--kind"]]
elif command == "clinvar-contract" and rest[:1] in (["pair-capture"], ["pair-check"], ["pair-complete"]):
    pairs = options(rest[1:], {"--vcf", "--producer", "--contract-source", "--aa", "--splice",
                               "--receipt", "--captured", "--config"})
    script = "clinvar_annotation_contract.py"
    forwarded = [rest[0], *sum((pair for pair in pairs if pair[0] != "--contract-source"), ())]
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
elif command == "mapped-artifact" and rest[:1] in (["run"], ["recover"]):
    head, producer = rest[1:], None
    if "--" in head:
        head, producer = head[:head.index("--")], head[head.index("--") + 1:]
    pairs = options(head, {"--artifact", "--mapping", "--input", "--cache-root", "--output-dir", "--assembly",
                           "--config", "--kind", "--helper-root", "--native-producer", "--expected-identity"})
    values = dict(pairs)
    if values.get("--native-producer") != "pm1-regions" or values.get("--kind") != "pm1":
        raise SystemExit("priva fixture: only the native PM1 producer has a Python equivalent")
    script = "mapped_artifact_publish.py"
    forwarded = [rest[0], *sum((pair for pair in pairs if pair[0] not in ("--helper-root", "--native-producer")), ()),
                 "--producer", os.path.join(SCRIPTS, "prepare_pm1_regions.py")]
    if producer is not None:
        built = dict(options(producer, {"--assembly", "--nm-to-enst", "--intolerant-domains", "--hcseeker",
                                        "--rmc", "--output"}))
        # The Python producer reads the DAS set and anchors from the config and checkout.
        forwarded += ["--", "--config", values["--config"], "--output", built["--output"],
                      "--assembly", built["--assembly"], "--nm_to_enst_tsv", built["--nm-to-enst"], "--skip_das"]
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
