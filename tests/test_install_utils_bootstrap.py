"""Installer startup tests with isolated fake managers; no packages/downloads.

Run with Python's standard library: python -m unittest discover -s tests
-p test_install_utils_bootstrap.py. Also discoverable by pytest.
"""

import json
import contextlib
import io
import os
from pathlib import Path
import shlex
import shutil
import subprocess
import sys
import tempfile
import textwrap
import unittest
from unittest import mock


ROOT = Path(__file__).resolve().parents[1]
INSTALLER = ROOT / "scripts" / "install_utils.sh"

FAKE_MANAGER = r'''
import json
import os
from pathlib import Path
import re
import shlex
import sys

root = Path(os.environ["PRIVA_TEST_ROOT"])
args = sys.argv[1:]
kind = Path(sys.argv[0]).name

def record(value):
    with (root / "calls.jsonl").open("a") as handle:
        handle.write(json.dumps(value) + "\n")

record([kind, *args])

def option(flag):
    return args[args.index(flag) + 1]

def create_environment(name, release=None):
    prefix = root / "envs" / name
    directory = prefix / "bin"
    directory.mkdir(parents=True, exist_ok=True)
    for tool in ("yq", "samtools", "vep", "vep_install", "perl"):
        if tool == os.environ.get("PRIVA_TEST_MISSING_TOOL"):
            continue
        destination = directory / tool
        destination.write_text(Path(os.environ["PRIVA_TEST_TOOL_TEMPLATE"]).read_text())
        destination.chmod(0o755)
    if not (directory / "python").exists() and os.environ.get("PRIVA_TEST_MISSING_TOOL") != "python":
        (directory / "python").write_text(Path(os.environ["PRIVA_TEST_PYTHON_TEMPLATE"]).read_text())
        (directory / "python").chmod(0o755)
    (prefix / "vep-release").write_text(release or os.environ.get("PRIVA_TEST_INITIAL_RELEASE", "113.2"))

if args == ["env", "list"]:
    if os.environ.get("PRIVA_TEST_FAIL_LIST"):
        sys.exit(41)
    print("Name Active Path")
    for directory in sorted((root / "envs").iterdir()):
        print(directory.name, "*", directory)
elif args[:2] == ["env", "create"]:
    if os.environ.get("PRIVA_TEST_FAIL_CREATE"):
        sys.exit(42)
    assert "--yes" in args
    specification = Path(option("--file")).read_text()
    record(["environment-spec", specification])
    requested = re.search(r"ensembl-vep=(\d+)\.\*", specification)
    create_environment(option("--name"), requested.group(1) + ".3" if requested else None)
elif args[:1] == ["install"]:
    assert "--yes" in args and "--freeze-installed" not in args
    pins = Path(option("--file")).read_text()
    record(["dependency-pins", pins])
    if os.environ.get("PRIVA_TEST_FAIL_INSTALL"):
        sys.exit(45)
    requested = next(re.fullmatch(r"ensembl-vep=(\d+)\.\*", arg).group(1)
                     for arg in args if re.fullmatch(r"ensembl-vep=(\d+)\.\*", arg))
    if not os.environ.get("PRIVA_TEST_IGNORE_INSTALL"):
        create_environment(option("--name"), requested + ".3")
elif args[:2] in (["shell", "activate"], ["shell.bash", "activate"]):
    if os.environ.get("PRIVA_TEST_FAIL_ACTIVATE"):
        sys.exit(43)
    if kind == "conda":
        assert args[:2] == ["shell.bash", "activate"]
        name = args[2]
    else:
        assert args[:2] == ["shell", "activate"]
        assert option("--shell") == "bash"
        name = option("--name")
    prefix = root / "envs" / name
    assert prefix.is_dir()
    print("export CONDA_PREFIX=" + shlex.quote(str(prefix)))
    print("export PATH=" + shlex.quote(str(prefix / "bin") + ":" + os.environ["PATH"]))
elif args[:1] == ["list"]:
    prefix = Path(option("--prefix")) if "--prefix" in args else root / "envs" / option("--name")
    assert prefix.is_dir()
    if "--json" in args:
        records = [
            {"name": "ensembl-vep", "version": (prefix / "vep-release").read_text(),
             "build_string": "hdfd78af_0", "channel": "bioconda"},
            {"name": "numpy", "version": "1.26.4", "build_string": "py312_0", "channel": "conda-forge"},
            {"name": "perl", "version": "5.32.1", "build": "7_hd590300_perl5", "channel": "conda-forge"},
            {"name": "pip-only", "version": "7.0", "build_string": "pypi_0", "channel": "pypi"},
        ]
        print(json.dumps({"packages": records} if os.environ.get("PRIVA_TEST_WRAPPED_LIST") else records))
    else:
        for name in ("ensembl-vep", "perl-bio-procedural", "perl-bioperl"):
            print(name, "1.0 test")
elif args[:4] == ["env", "config", "vars", "set"]:
    assert option("--prefix") == os.environ["CONDA_PREFIX"]
    if os.environ.get("PRIVA_TEST_FAIL_VARS"):
        sys.exit(44)
else:
    raise SystemExit("Unexpected fake manager call: " + repr(args))
'''


FAKE_PYTHON = r'''
import json
import os
from pathlib import Path
import sys

# Bootstrap fixtures do not depend on the host's scientific packages.
if sys.argv[1:] == ["-I", "-c", "import diptest, pyarrow"]:
    with (Path(os.environ["PRIVA_TEST_ROOT"]) / "calls.jsonl").open("a") as handle:
        handle.write(json.dumps(["python-import-check", *sys.argv[1:]]) + "\n")
    failure = os.environ.get("PRIVA_TEST_IMPORT_FAILURE")
    if failure:
        print(failure, file=sys.stderr)
        sys.exit(71)
    sys.exit(0)
os.execv(sys.executable, [sys.executable, *sys.argv[1:]])
'''

FAKE_TOOL = r'''
import json
import os
from pathlib import Path
import re
import sys

root = Path(os.environ["PRIVA_TEST_ROOT"])
tool = Path(sys.argv[0]).name
args = sys.argv[1:]
prefix = Path(os.environ["CONDA_PREFIX"])
release = (prefix / "vep-release").read_text()
major = release.split(".")[0]
with (root / "calls.jsonl").open("a") as handle:
    handle.write(json.dumps([tool, *args]) + "\n")

def option(flag):
    return args[args.index(flag) + 1]

if tool == "yq":
    assert prefix / "bin" == Path(sys.argv[0]).parent
    text = Path(sys.argv[2]).read_text()
    try:
        config = json.loads(text)
    except ValueError:
        config = {}
        for line in text.splitlines():
            match = re.match(r"^([a-z_]+):\s*(.*?)\s*(?:#.*)?$", line)
            if match:
                value = match.group(2)
                try:
                    value = json.loads(value)
                except ValueError:
                    value = value.strip("'")
                config[match.group(1)] = value
    print(json.dumps(config if sys.argv[1] == "." else config.get(sys.argv[1].lstrip("."))))
    if os.environ.get("PRIVA_TEST_YQ_FAIL_WITH_VALUE"):
        sys.exit(2)
elif tool == "vep":
    if os.environ.get("PRIVA_TEST_FAIL_VEP"):
        sys.exit(39)
    print("ensembl-vep", os.environ.get("PRIVA_TEST_REPORTED_RELEASE", release))
elif tool == "perl":
    if "-MBio::EnsEMBL::ApiVersion" in args:
        print(os.environ.get("PRIVA_TEST_API_RELEASE", major), end="")
    elif "-MBio::EnsEMBL::VEP::Utils=get_version_data" in args:
        print(os.environ.get("PRIVA_TEST_HELPER_RELEASE", major), end="")
    else:
        print(str(prefix / "lib/site_perl"))
elif tool == "vep_install":
    if os.environ.get("PRIVA_TEST_INSTALLER_FAIL"):
        sys.exit(37)
    automatic = option("--AUTO")
    selected_release = option("--CACHE_VERSION")
    if "c" in automatic and not os.environ.get("PRIVA_TEST_SKIP_CACHE"):
        cache = Path(option("--CACHEDIR")) / "homo_sapiens_merged" / (selected_release + "_" + option("--ASSEMBLY"))
        (cache / "1").mkdir(parents=True)
        (cache / "1/all_vars.gz").write_text("fixture cache, not genomic data\n")
        if not os.environ.get("PRIVA_TEST_OMIT_CACHE_INDEX"):
            (cache / ("1/all_vars.gz." + os.environ.get("PRIVA_TEST_CACHE_INDEX_SUFFIX", "csi"))).write_text("fixture index\n")
        (cache / "chr_synonyms.txt").write_text("fixture synonyms\n")
        assembly = os.environ.get("PRIVA_TEST_CACHE_ASSEMBLY", option("--ASSEMBLY"))
        species = os.environ.get("PRIVA_TEST_CACHE_SPECIES", "homo_sapiens")
        (cache / "info.txt").write_text("assembly\t" + assembly + "\nspecies\t" + species + "\nvar_type\ttabix\n")
    if "p" in automatic:
        directory = Path(option("--PLUGINSDIR"))
        (directory / "plugin_config.txt").write_text("official fixture release " + selected_release + "\n")
        for name in option("--PLUGINS").split(","):
            if name != os.environ.get("PRIVA_TEST_OMIT_PLUGIN"):
                if name == "UTRAnnotator" and selected_release == "116":
                    (directory / (name + ".pm")).write_bytes(Path(os.environ["PRIVA_TEST_UTR116_ORIGINAL"]).read_bytes())
                else:
                    (directory / (name + ".pm")).write_text("package " + name + ";\n# official fixture release " + selected_release + "\n1;\n")
'''


class InstallerBootstrapTests(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory(prefix="priva_bootstrap_", dir=os.environ.get("PRIVA_TEST_TMPDIR", os.environ.get("TMPDIR")))
        self.addCleanup(self.temporary.cleanup)
        self.directory = Path(self.temporary.name)
        self.bin = self.directory / "bin"
        self.bin.mkdir()
        (self.directory / "envs").mkdir()
        self.bash = shutil.which("bash")
        for command in ("bash", "sed", "awk", "dirname", "realpath", "date", "basename", "grep", "head", "mktemp", "rm", "mkdir", "getopt"):
            executable = shutil.which(command)
            self.assertIsNotNone(executable, command)
            (self.bin / command).symlink_to(executable)
        self.tool_template = self.directory / "tool_template.py"
        self.write_executable(self.tool_template, FAKE_TOOL)
        self.python_template = self.directory / "python_template.py"
        self.write_executable(self.python_template, FAKE_PYTHON)
        self.config = self.directory / "config.yaml"
        (self.directory / "selected_cache").mkdir()
        self.config.write_text("".join(f"{key}: {json.dumps(value)}\n" for key, value in {
            "vep_installed_version": "113",
            "cache_root": str(self.directory / "selected_cache"), "tmp_dir": str(self.directory),
            "assembly": "hg38", "ref_genome": str(self.directory / "reference.fa"),
            "vep_cache_dir": str(self.directory), "vep_plugins_dir": str(self.directory),
            "vep_plugins_cachedir": str(self.directory),
            "conda_env_yaml": "/old/machine/ignored_environment.yml",
        }.items()))
        self.env = {key: value for key, value in os.environ.items()
                    if not key.startswith(("CONDA", "MAMBA", "PRIVA_", "BASH_FUNC_"))
                    and key not in ("BASH_ENV", "ENV")}
        self.env.update({
            "PATH": str(self.bin), "PRIVA_TEST_ROOT": str(self.directory),
            "PRIVA_TEST_TOOL_TEMPLATE": str(self.tool_template),
            "PRIVA_TEST_PYTHON_TEMPLATE": str(self.python_template),
            "PRIVA_TEST_UTR116_ORIGINAL": str(ROOT / "tests/fixtures/UTRAnnotator116.original.pm"),
            "PRIVA_SMOKE_SKIP_VEP": "1", "TMPDIR": str(self.directory),
        })

    def write_executable(self, path, body):
        path.write_text("#!" + sys.executable + "\n" + textwrap.dedent(body))
        path.chmod(0o755)

    def manager(self, kind, outside_path=False):
        parent = self.directory if outside_path else self.bin
        path = parent / kind
        self.write_executable(path, FAKE_MANAGER)
        if outside_path:
            self.env["CONDA_EXE" if kind == "conda" else "MAMBA_EXE"] = str(path)
        return path

    def existing_environment(self, name="priva_acmg", missing_tool=None):
        directory = self.directory / "envs" / name / "bin"
        directory.mkdir(parents=True)
        for tool in ("yq", "samtools", "vep", "vep_install", "perl"):
            if tool != missing_tool:
                self.write_executable(directory / tool, FAKE_TOOL)
        if missing_tool != "python":
            self.write_executable(directory / "python", FAKE_PYTHON)
        (directory.parent / "vep-release").write_text(self.env.get("PRIVA_TEST_INITIAL_RELEASE", "113.2"))

    def calls(self):
        path = self.directory / "calls.jsonl"
        return [json.loads(line) for line in path.read_text().splitlines()] if path.exists() else []

    def run_shell(self, body, *args):
        script = 'source "$1" >/dev/null || exit; shift\n' + body
        return subprocess.run(
            [self.bash, "--noprofile", "--norc", "-c", script, "bash", str(INSTALLER), *map(str, args)],
            env=self.env, capture_output=True, text=True, timeout=15,
        )

    def run_main(self, *args):
        # Stop at the first real resource operation, after startup/config reads.
        return self.run_shell('''
            # Resource preparation is tested separately with actual sources.
            # These environment-manager fixtures deliberately have no reference
            # or external caches and must not invoke source adoption.
            priva_prepare_install_config() { return 0; }
            bundled_snapshot_install() { return 0; }
            Nextflow_runtime_install() { return 0; }  # Driver acquisition has separate fixtures.
            PriVA_rust_install() { return 0; }  # Rust deployment has separate fixtures.
            reference_genome_install() { printf 'REFERENCE_STAGE_REACHED\n'; return 73; }
            main_install "$@"
        ''', *args)

    def test_bundled_environment_bootstraps_without_yq_for_each_manager(self):
        for kind in ("micromamba", "mamba", "conda"):
            with self.subTest(manager=kind):
                manager = self.manager(kind)
                result = self.run_main(self.config)
                self.assertIn("REFERENCE_STAGE_REACHED", result.stdout, result.stderr)
                calls = self.calls()
                create = next(call for call in calls if call[1:3] == ["env", "create"])
                self.assertNotEqual(create[create.index("--file") + 1], str(ROOT / "acmg_conda_rough.yml"))
                self.assertFalse(Path(create[create.index("--file") + 1]).exists())
                specification = next(call[1] for call in calls if call[0] == "environment-spec")
                self.assertIn("ensembl-vep=113.*", specification)
                self.assertIn("  - diptest\n", specification)
                self.assertIn("  - pyarrow\n", specification)
                self.assertEqual(create[create.index("--name") + 1], "priva_acmg")
                first_yq = next(index for index, call in enumerate(calls) if call[0] == "yq")
                activation = next(index for index, call in enumerate(calls) if call[1:3] in
                                  (["shell", "activate"], ["shell.bash", "activate"]))
                self.assertLess(activation, first_yq)
                manager.unlink()
                shutil.rmtree(self.directory / "envs" / "priva_acmg")
                (self.directory / "calls.jsonl").unlink()

    def test_environment_only_entry_point_needs_no_file_argument(self):
        self.manager("micromamba")
        result = self.run_shell('conda_install_vep')
        self.assertEqual(result.returncode, 0, result.stderr)
        create = next(call for call in self.calls() if call[1:3] == ["env", "create"])
        self.assertEqual(create[create.index("--file") + 1], str(ROOT / "acmg_conda_rough.yml"))
        self.assertEqual(create[create.index("--name") + 1], "priva_acmg")
        self.assertFalse(any(call[0] == "yq" for call in self.calls()))

    def test_environment_name_stays_fixed_without_parsing_yaml(self):
        self.manager("conda")
        env_file = self.directory / "acmg_conda_rough.yml"
        env_file.write_text("name: another_name\ndependencies: [yq]\n")
        result = self.run_shell('BASE_DIR=$1; conda_install_vep', self.directory)
        self.assertEqual(result.returncode, 0, result.stderr)
        create = next(call for call in self.calls() if call[1:3] == ["env", "create"])
        self.assertEqual(create[create.index("--file") + 1], str(env_file))
        self.assertEqual(create[create.index("--name") + 1], "priva_acmg")

    def test_existing_environment_reused_and_wrapper_found_outside_path(self):
        self.manager("micromamba", outside_path=True)
        self.existing_environment()
        result = self.run_main(self.config)
        self.assertIn("REFERENCE_STAGE_REACHED", result.stdout, result.stderr)
        self.assertFalse(any(call[1:3] == ["env", "create"] for call in self.calls()))

    def test_similar_environment_name_does_not_skip_creation(self):
        self.manager("mamba")
        self.existing_environment("priva_acmg_backup")
        result = self.run_main(self.config)
        self.assertIn("REFERENCE_STAGE_REACHED", result.stdout, result.stderr)
        self.assertTrue(any(call[1:3] == ["env", "create"] for call in self.calls()))

    def test_conda_backed_mamba_uses_sibling_conda_activation(self):
        self.manager("mamba")
        self.manager("conda")
        result = self.run_main(self.config)
        self.assertIn("REFERENCE_STAGE_REACHED", result.stdout, result.stderr)
        self.assertTrue(any(call[:3] == ["mamba", "env", "create"] for call in self.calls()))
        self.assertTrue(any(call[:3] == ["conda", "shell.bash", "activate"] for call in self.calls()))

    def test_help_and_invalid_arguments_do_not_touch_a_manager(self):
        self.manager("micromamba")
        direct_help = subprocess.run(
            [self.bash, str(INSTALLER), "main_install", "--help"],
            env=self.env, capture_output=True, text=True, timeout=15,
        )
        self.assertEqual(direct_help.returncode, 0, direct_help.stderr)
        self.assertIn("priva_acmg", direct_help.stdout)
        self.assertNotIn("--env-file", direct_help.stdout)
        help_result = self.run_main("--help")
        self.assertEqual(help_result.returncode, 0, help_result.stderr)
        self.assertIn("priva_acmg", help_result.stdout)
        self.assertNotIn("--env-file", help_result.stdout)
        for args in ((), (self.config, "--env-file", ROOT / "acmg_conda.yml"),
                     (self.config, "--env-file=" + str(ROOT / "acmg_conda.yml")),
                     ("--env-file", ROOT / "acmg_conda.yml", self.config),
                     (self.config, "--unknown"), (self.config, self.config),
                     (self.directory / "missing_config",), ("--help", self.config)):
            with self.subTest(args=args):
                self.assertNotEqual(self.run_main(*args).returncode, 0)
        self.assertNotEqual(self.run_shell('conda_install_vep "$@"', ROOT / "acmg_conda.yml").returncode, 0)
        self.assertEqual(self.calls(), [])

    def test_missing_manager_fails_before_reading_config(self):
        result = self.run_main(self.config)
        self.assertNotEqual(result.returncode, 0)
        self.assertIn("Install micromamba, mamba, or conda", result.stderr)
        self.assertEqual(self.calls(), [])

    def test_missing_bundled_file_fails_before_manager_operations(self):
        self.manager("micromamba")
        result = self.run_shell('BASE_DIR=$1; conda_install_vep', self.directory)
        self.assertNotEqual(result.returncode, 0)
        self.assertIn("Environment file not readable", result.stderr)
        self.assertEqual(self.calls(), [])

    def test_manager_failures_stop_before_yq_or_resource_work(self):
        self.manager("micromamba")
        for variable in ("PRIVA_TEST_FAIL_LIST", "PRIVA_TEST_FAIL_CREATE", "PRIVA_TEST_FAIL_ACTIVATE"):
            with self.subTest(variable=variable):
                self.env[variable] = "1"
                result = self.run_main(self.config)
                self.assertNotEqual(result.returncode, 0)
                self.assertNotIn("REFERENCE_STAGE_REACHED", result.stdout)
                self.assertFalse(any(call[0] == "yq" for call in self.calls()))
                self.env.pop(variable)

    def test_existing_environment_missing_yq_is_not_accepted(self):
        self.manager("micromamba")
        self.existing_environment(missing_tool="yq")
        result = self.run_main(self.config)
        self.assertNotEqual(result.returncode, 0)
        self.assertIn("missing its own yq", result.stderr)
        self.assertFalse(any(call[1:3] == ["env", "create"] for call in self.calls()))

    def test_required_python_imports_are_checked_before_resource_work(self):
        self.manager("micromamba")
        self.existing_environment()
        for failure in ("ModuleNotFoundError: diptest", "ModuleNotFoundError: pyarrow",
                        "ImportError: incompatible binary ABI"):
            with self.subTest(failure=failure):
                self.env["PRIVA_TEST_IMPORT_FAILURE"] = failure
                result = self.run_main(self.config)
                self.assertNotEqual(result.returncode, 0)
                self.assertNotIn("REFERENCE_STAGE_REACHED", result.stdout)
                self.assertIn(failure, result.stderr)
                self.assertIn("cannot import required diptest/pyarrow", result.stderr)
                self.assertFalse(any(call[1:3] == ["env", "create"] or call[1:2] == ["install"]
                                     for call in self.calls()))
        self.env.pop("PRIVA_TEST_IMPORT_FAILURE")
        result = self.run_main(self.config)
        self.assertIn("REFERENCE_STAGE_REACHED", result.stdout, result.stderr)
        self.assertTrue(any(call[0] == "python-import-check" for call in self.calls()))

    def test_existing_environment_missing_own_python_stops_before_resource_work(self):
        self.manager("micromamba")
        self.existing_environment(missing_tool="python")
        result = self.run_main(self.config)
        self.assertNotEqual(result.returncode, 0)
        self.assertIn("missing its own python", result.stderr)
        self.assertNotIn("REFERENCE_STAGE_REACHED", result.stdout)

    def test_saved_variables_use_selected_manager_and_apply_immediately(self):
        self.manager("micromamba")
        result = self.run_shell('''
            conda_install_vep || exit
            priva_env_vars_set 'PERL5LIB=/example/perl path' || exit
            printf 'PERL5LIB=%s\n' "$PERL5LIB"
        ''')
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertIn("PERL5LIB=/example/perl path", result.stdout)
        self.assertTrue(any(call[:5] == ["micromamba", "env", "config", "vars", "set"] for call in self.calls()))

    def select_release(self, release):
        self.config.write_text(self.config.read_text().replace(
            'vep_installed_version: "113"', f'vep_installed_version: "{release}"'))

    def cache_paths(self):
        cache = self.directory / "cache"
        plugins = self.directory / "plugins"
        cache.mkdir(exist_ok=True)
        plugins.mkdir(exist_ok=True)
        self.config.write_text(self.config.read_text().replace(
            "vep_cache_dir: " + json.dumps(str(self.directory)),
            "vep_cache_dir: " + json.dumps(str(cache))).replace(
            "vep_plugins_dir: " + json.dumps(str(self.directory)),
            "vep_plugins_dir: " + json.dumps(str(plugins))).replace(
            "cache_root: " + json.dumps(str(self.directory / "selected_cache")),
            "cache_root: " + json.dumps(str(self.directory))))
        if "cache_folders:" not in self.config.read_text():
            with self.config.open("a") as handle:
                handle.write('cache_folders: {"vep_plugin": "plugins"}\n')
        return cache, plugins

    def run_cache(self, assembly="hg19", release="113"):
        cache, plugins = self.cache_paths()
        return self.run_shell(r'''
            conda_install_vep "$1" || exit
            vep_cache_api_install --VEP_VERSION "$2" --VEP_ASSEMBLY "$3" \
                --VEP_CACHEDIR "$4" --VEP_PLUGINSDIR "$5"
        ''', self.config, release, assembly, cache, plugins)

    def test_vep_failure_is_not_masked_by_environment_setup(self):
        self.manager("micromamba")
        self.env["PRIVA_TEST_INSTALLER_FAIL"] = "1"
        result = self.run_cache()
        self.assertNotEqual(result.returncode, 0)
        self.assertIn("Failed to install the VEP cache/plugins", result.stderr)
        self.assertFalse(any(call[1:5] == ["env", "config", "vars", "set"] for call in self.calls()))
        self.assertFalse(list(self.directory.glob("plugins.install-*")))

    def test_installer_retry_recovers_utr_receipt_failure_before_regular_plugin_install(self):
        import importlib.util
        spec = importlib.util.spec_from_file_location("utr_retry_fixture", ROOT / "scripts/vep_plugins_release.py")
        helper = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(helper)
        self.select_release("116")
        self.manager("micromamba")
        first = self.run_cache(release="116")
        self.assertEqual(first.returncode, 0, first.stderr)
        plugins = self.directory / "plugins"
        target, receipt_file = plugins / "UTRAnnotator.pm", plugins / helper.RECEIPT
        receipt = helper.load_receipt(plugins)
        original = Path(receipt["local_patches"].pop("UTRAnnotator.pm")["original_file"])
        receipt.pop("local_patches")
        target.write_bytes(original.read_bytes())
        helper.write_receipt(receipt_file, receipt)
        write = helper.write_receipt
        def fail_receipt(path, value):
            if path == receipt_file:
                raise OSError("simulated live receipt failure")
            return write(path, value)
        with mock.patch.object(helper, "write_receipt", side_effect=fail_receipt):
            with self.assertRaises(OSError):
                helper.ensure_utr_reporting(plugins, "116")
        self.assertFalse(helper.matches(plugins, helper.load_receipt(plugins), "116", ["UTRAnnotator"]))
        before = len([call for call in self.calls() if call[0] == "vep_install"])
        retried = self.run_cache(release="116")
        self.assertEqual(retried.returncode, 0, retried.stderr)
        self.assertEqual(len([call for call in self.calls() if call[0] == "vep_install"]), before)
        self.assertFalse((plugins / ".priva-utr-reporting.pending.json").exists())
        self.assertEqual(helper.digest(target), helper.UTR_INSTALLED_SHA)
        self.assertEqual(helper.digest(original), helper.UTR_UPSTREAM_SHA)

    def test_malformed_utr_recovery_stops_installer_before_live_replacement(self):
        self.select_release("116")
        self.manager("micromamba")
        first = self.run_cache(release="116")
        self.assertEqual(first.returncode, 0, first.stderr)
        plugins = self.directory / "plugins"
        pending = plugins / ".priva-utr-reporting.pending.json"
        pending.write_text("null\n")
        watched = (pending, plugins / "UTRAnnotator.pm", plugins / ".priva-vep-plugins.json")
        before = [path.read_bytes() for path in watched]
        calls = len([call for call in self.calls() if call[0] == "vep_install"])
        retried = self.run_cache(release="116")
        self.assertNotEqual(retried.returncode, 0)
        self.assertIn("recovery receipt", retried.stderr)
        self.assertEqual([path.read_bytes() for path in watched], before)
        self.assertEqual(len([call for call in self.calls() if call[0] == "vep_install"]), calls)

    def test_vep_zero_exit_without_cache_is_not_success(self):
        self.manager("micromamba")
        self.env["PRIVA_TEST_SKIP_CACHE"] = "1"
        result = self.run_cache()
        self.assertNotEqual(result.returncode, 0)
        install = next(call for call in self.calls() if call[0] == "vep_install")
        self.assertEqual(install[install.index("--AUTO") + 1], "cp")
        self.assertEqual(install[install.index("--CACHE_VERSION") + 1], "113")
        self.assertIn("without required cache files", result.stderr)
        self.assertFalse((self.directory / "plugins/.priva-vep-plugins.json").exists())
        self.assertFalse(any(call[1:5] == ["env", "config", "vars", "set"] for call in self.calls()))

    def test_fresh_requested_release_is_pinned_for_each_manager(self):
        self.select_release("116")
        for kind in ("micromamba", "mamba", "conda"):
            with self.subTest(manager=kind):
                manager = self.manager(kind)
                result = self.run_main(self.config)
                self.assertIn("REFERENCE_STAGE_REACHED", result.stdout, result.stderr)
                specification = next(call[1] for call in self.calls() if call[0] == "environment-spec")
                self.assertIn("ensembl-vep=116.*", specification)
                self.assertEqual((self.directory / "envs/priva_acmg/vep-release").read_text(), "116.3")
                self.assertIn('vep_installed_version: "116"', self.config.read_text())
                manager.unlink()
                shutil.rmtree(self.directory / "envs/priva_acmg")
                (self.directory / "calls.jsonl").unlink()

    def test_existing_release_switch_pins_other_packages_and_supports_wrapped_json(self):
        self.manager("micromamba")
        self.existing_environment()
        self.select_release("116")
        self.env["PRIVA_TEST_WRAPPED_LIST"] = "1"
        result = self.run_main(self.config)
        self.assertIn("REFERENCE_STAGE_REACHED", result.stdout, result.stderr)
        install = next(call for call in self.calls() if call[1:2] == ["install"])
        self.assertIn("ensembl-vep=116.*", install)
        self.assertNotIn("--freeze-installed", install)
        self.assertFalse(Path(install[install.index("--file") + 1]).exists())
        pins = next(call[1] for call in self.calls() if call[0] == "dependency-pins")
        self.assertEqual(set(pins.splitlines()), {"numpy=1.26.4=py312_0", "perl=5.32.1=7_hd590300_perl5"})
        package_list = next(call for call in self.calls() if "--json" in call)
        self.assertEqual(package_list[package_list.index("--prefix") + 1],
                         str(self.directory / "envs/priva_acmg"))
        self.assertIn('vep_installed_version: "116"', self.config.read_text())

    def test_existing_environment_can_deliberately_downgrade(self):
        self.manager("conda")
        self.env["PRIVA_TEST_INITIAL_RELEASE"] = "116.3"
        self.existing_environment()
        result = self.run_main(self.config)
        self.assertIn("REFERENCE_STAGE_REACHED", result.stdout, result.stderr)
        self.assertTrue(any("ensembl-vep=113.*" in call for call in self.calls()))
        self.assertEqual((self.directory / "envs/priva_acmg/vep-release").read_text(), "113.3")

    def test_same_major_existing_environment_is_not_updated(self):
        self.manager("micromamba")
        self.existing_environment()
        result = self.run_main(self.config)
        self.assertIn("REFERENCE_STAGE_REACHED", result.stdout, result.stderr)
        self.assertFalse(any(call[1:2] == ["install"] for call in self.calls()))

    def test_invalid_release_fails_before_manager_mutation(self):
        self.manager("micromamba")
        invalid = ("", "vep_installed_version: 116.3\n", "vep_installed_version: latest\n",
                   "vep_installed_version: $(touch SHOULD_NOT_EXIST)\n", "vep_installed_version: 0116\n",
                   "vep_installed_version: true\n", 'vep_installed_version: "113"\nvep_installed_version: 116\n',
                   'vep_installed_version: "113"\n"vep_installed_version": 116\n',
                   'vep_installed_version: 113\n---\nvep_installed_version: 116\n',
                   'vep_installed_version: 113\n--- {vep_installed_version: 116}\n',
                   '{"vep_installed_version": 116}\n')
        for text in invalid:
            with self.subTest(config=text):
                self.config.write_text(text)
                result = self.run_main(self.config)
                self.assertNotEqual(result.returncode, 0)
                self.assertNotIn("REFERENCE_STAGE_REACHED", result.stdout)
                self.assertEqual(self.calls(), [])

    def test_numeric_release_scalar_accepts_quotes_and_comments(self):
        self.manager("micromamba")
        for scalar in ("116", "'116'", '"116"'):
            with self.subTest(scalar=scalar):
                self.config.write_text("vep_installed_version: " + scalar + " # requested\n")
                result = self.run_shell('conda_install_vep "$1"', self.config)
                self.assertEqual(result.returncode, 0, result.stderr)
                shutil.rmtree(self.directory / "envs/priva_acmg")
                (self.directory / "calls.jsonl").unlink()

    def test_unavailable_or_conflicting_release_preserves_requested_config(self):
        self.manager("micromamba")
        self.existing_environment()
        self.select_release("999")
        self.env["PRIVA_TEST_FAIL_INSTALL"] = "1"
        before = self.config.read_bytes()
        result = self.run_main(self.config)
        self.assertNotEqual(result.returncode, 0)
        self.assertIn("unavailable or incompatible", result.stderr)
        self.assertNotIn("REFERENCE_STAGE_REACHED", result.stdout)
        self.assertEqual(self.config.read_bytes(), before)
        self.assertEqual((self.directory / "envs/priva_acmg/vep-release").read_text(), "113.2")
        self.assertFalse(list(self.directory.glob("priva-vep-pins.*")))

    def test_solver_success_with_wrong_release_fails_verification(self):
        self.manager("micromamba")
        self.existing_environment()
        self.select_release("116")
        self.env["PRIVA_TEST_IGNORE_INSTALL"] = "1"
        result = self.run_main(self.config)
        self.assertNotEqual(result.returncode, 0)
        self.assertIn("Requested VEP 116; executable=113", result.stderr)
        self.assertNotIn("REFERENCE_STAGE_REACHED", result.stdout)

    def test_yaml_parse_failure_cannot_be_hidden_before_existing_environment_switch(self):
        self.manager("micromamba")
        self.existing_environment()
        self.select_release("116")
        self.env["PRIVA_TEST_YQ_FAIL_WITH_VALUE"] = "1"
        result = self.run_main(self.config)
        self.assertNotEqual(result.returncode, 0)
        self.assertIn("Could not parse the installer configuration", result.stderr)
        self.assertFalse(any(call[1:2] == ["install"] for call in self.calls()))
        self.assertEqual((self.directory / "envs/priva_acmg/vep-release").read_text(), "113.2")

    def test_main_forwards_requested_release_without_rewriting_config(self):
        self.manager("micromamba")
        self.select_release("116")
        self.env["PRIVA_SMOKE_SKIP_VEP"] = "0"
        (self.directory / "reference.fa").write_text(">1\nA\n")
        before = self.config.read_bytes()
        result = self.run_shell(r'''
            # Release-forwarding boundary fixture, not a manual-cache install.
            priva_manual_cache_inputs() { return 0; }
            priva_install_prerequisites() { return 0; }
            nm_to_enst_mapping_install() { return 0; }
            priva_configure_cache_paths() { [[ ${2:-} != --check ]] || printf '%s' '{}'; return 0; }
            priva_configure_portable_paths() { return 0; }
            Nextflow_runtime_install() { return 0; }  # Release-forwarding fixture only.
            PriVA_rust_install() { return 0; }  # Release-forwarding fixture only.
            bundled_snapshot_install() { return 0; }  # Release-forwarding fixture only.
            reference_genome_install() { return 0; }
            update_yaml() { printf 'UNEXPECTED_CONFIG_UPDATE\n'; return 1; }
            vep_cache_api_install() { printf 'CACHE_ARGS=%s\n' "$*"; return 73; }
            main_install "$1"
        ''', self.config)
        self.assertNotEqual(result.returncode, 0)
        self.assertIn("--VEP_VERSION 116", result.stdout)
        self.assertNotIn("UNEXPECTED_CONFIG_UPDATE", result.stdout)
        self.assertEqual(self.config.read_bytes(), before)

    def test_api_or_installer_mismatch_stops_before_resources(self):
        self.manager("micromamba")
        self.existing_environment()
        for variable in ("PRIVA_TEST_API_RELEASE", "PRIVA_TEST_HELPER_RELEASE"):
            with self.subTest(variable=variable):
                self.env[variable] = "116"
                result = self.run_main(self.config)
                self.assertNotEqual(result.returncode, 0)
                self.assertIn("Requested VEP 113", result.stderr)
                self.assertNotIn("REFERENCE_STAGE_REACHED", result.stdout)
                self.env.pop(variable)

    def test_caches_use_requested_major_and_both_assembly_aliases(self):
        self.manager("micromamba")
        self.select_release("116")
        for assembly, expected in (("hg19", "GRCh37"), ("hg38", "GRCh38")):
            with self.subTest(assembly=assembly):
                result = self.run_cache(assembly, "116")
                self.assertEqual(result.returncode, 0, result.stderr)
                cache = self.directory / ("cache/homo_sapiens_merged/116_" + expected)
                self.assertTrue((cache / "1/all_vars.gz").is_file())
                receipt = json.loads((self.directory / "plugins/.priva-vep-plugins.json").read_text())
                self.assertEqual(receipt["release"], "116")
                self.assertEqual(receipt["source"], "https://raw.githubusercontent.com/Ensembl/VEP_plugins/release/116")
                self.assertEqual(len(receipt["plugins"]), 12)
                self.assertIn('vep_installed_version: "116"', self.config.read_text())

    def test_existing_incomplete_cache_fails_without_deleting_or_downloading(self):
        self.manager("micromamba")
        cache, _ = self.cache_paths()
        directory = cache / "homo_sapiens_merged/113_GRCh37"
        directory.mkdir(parents=True)
        sentinel = directory / "partial"
        sentinel.write_text("retained incomplete manual cache")
        result = self.run_cache()
        self.assertNotEqual(result.returncode, 0)
        self.assertIn("incomplete or mismatched", result.stderr)
        self.assertEqual(sentinel.read_text(), "retained incomplete manual cache")
        self.assertFalse(any(call[0] == "vep_install" for call in self.calls()))

    def test_cache_conversion_warning_without_index_is_not_success(self):
        self.manager("micromamba")
        self.env["PRIVA_TEST_OMIT_CACHE_INDEX"] = "1"
        result = self.run_cache()
        self.assertNotEqual(result.returncode, 0)
        self.assertIn("without required cache files", result.stderr)
        self.assertFalse((self.directory / "plugins/.priva-vep-plugins.json").exists())

    def test_converted_cache_accepts_tbi_as_well_as_csi_indexes(self):
        self.manager("micromamba")
        self.env["PRIVA_TEST_CACHE_INDEX_SUFFIX"] = "tbi"
        result = self.run_cache()
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertTrue((self.directory / "cache/homo_sapiens_merged/113_GRCh37/1/all_vars.gz.tbi").exists())

    def test_downloaded_wrong_assembly_or_species_cannot_publish_plugins(self):
        self.manager("micromamba")
        for variable, value in (("PRIVA_TEST_CACHE_ASSEMBLY", "GRCh38"),
                                ("PRIVA_TEST_CACHE_SPECIES", "mus_musculus")):
            with self.subTest(variable=variable):
                self.env[variable] = value
                result = self.run_cache()
                self.assertNotEqual(result.returncode, 0)
                self.assertIn("matching assembly/species", result.stderr)
                self.assertFalse((self.directory / "plugins/.priva-vep-plugins.json").exists())
                shutil.rmtree(self.directory / "cache")
                self.env.pop(variable)

    def test_missing_plugin_never_uses_stale_existing_file(self):
        self.manager("micromamba")
        _, plugins = self.cache_paths()
        stale = plugins / "SpliceAI.pm"
        stale.write_text("package SpliceAI;\n# old customized plugin\n1;\n")
        self.env["PRIVA_TEST_OMIT_PLUGIN"] = "SpliceAI"
        result = self.run_cache()
        self.assertNotEqual(result.returncode, 0)
        self.assertIn("Missing or empty regular plugin", result.stderr)
        self.assertIn("old customized", stale.read_text())
        self.assertFalse((plugins / ".priva-vep-plugins.json").exists())

    def test_same_major_customized_plugin_is_preserved(self):
        self.manager("micromamba")
        _, plugins = self.cache_paths()
        custom = plugins / "SpliceAI.pm"
        custom.write_text("package SpliceAI;\n# custom\n1;\n")
        result = self.run_cache()
        self.assertNotEqual(result.returncode, 0)
        self.assertIn("differs from official release 113", result.stderr)
        self.assertIn("# custom", custom.read_text())
        self.assertFalse((plugins / ".priva-vep-plugins.json").exists())
        self.assertFalse((plugins / "SpliceVault.pm").exists())

    def test_verified_receipt_reuses_plugins_and_detects_changed_bytes(self):
        self.manager("micromamba")
        result = self.run_cache()
        self.assertEqual(result.returncode, 0, result.stderr)
        count = sum(call[0] == "vep_install" for call in self.calls())
        result = self.run_cache()
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertEqual(sum(call[0] == "vep_install" for call in self.calls()), count)
        plugin = self.directory / "plugins/SpliceAI.pm"
        plugin.write_text("package SpliceAI;\n# locally changed\n1;\n")
        result = self.run_cache()
        self.assertNotEqual(result.returncode, 0)
        self.assertIn("locally changed", plugin.read_text())
        self.assertEqual(sum(call[0] == "vep_install" for call in self.calls()), count + 1)

    def test_explicit_release_switch_refreshes_selected_plugins_preserving_loftee(self):
        self.manager("micromamba")
        self.existing_environment()
        _, plugins = self.cache_paths()
        (plugins / "SpliceAI.pm").write_text("package SpliceAI;\n# old release\n1;\n")
        loftee = plugins / "LoFtee_repos/loftee-hg19"
        loftee.mkdir(parents=True)
        sentinel = loftee / "LoF.pm"
        sentinel.write_text("external LOFTEE retained")
        self.select_release("116")
        result = self.run_cache(release="116")
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertIn("official fixture release 116", (plugins / "SpliceAI.pm").read_text())
        self.assertEqual(sentinel.read_text(), "external LOFTEE retained")

    def test_first_unrecorded_upgrade_retries_after_cache_failure(self):
        self.manager("micromamba")
        self.existing_environment()
        _, plugins = self.cache_paths()
        old = plugins / "SpliceAI.pm"
        old.write_text("package SpliceAI;\n# unrecorded old113\n1;\n")
        self.select_release("116")
        self.env["PRIVA_TEST_SKIP_CACHE"] = "1"
        result = self.run_cache(release="116")
        self.assertNotEqual(result.returncode, 0)
        self.assertIn("unrecorded old113", old.read_text())
        transition = self.directory / "envs/priva_acmg/.priva-vep-release-transition.json"
        self.assertTrue(transition.is_file())
        self.assertEqual(json.loads(transition.read_text())["from_release"], "113")
        self.assertEqual((self.directory / "envs/priva_acmg/vep-release").read_text(), "116.3")
        self.env.pop("PRIVA_TEST_SKIP_CACHE")
        result = self.run_cache(release="116")
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertIn("official fixture release 116", old.read_text())
        self.assertFalse(transition.exists())

    def test_unrecorded_upgrade_retry_preserves_new_custom_bytes(self):
        self.manager("micromamba")
        self.existing_environment()
        _, plugins = self.cache_paths()
        old = plugins / "SpliceAI.pm"
        old.write_text("package SpliceAI;\n# before upgrade\n1;\n")
        self.select_release("116")
        self.env["PRIVA_TEST_SKIP_CACHE"] = "1"
        self.assertNotEqual(self.run_cache(release="116").returncode, 0)
        old.write_text("package SpliceAI;\n# changed after failure\n1;\n")
        self.env.pop("PRIVA_TEST_SKIP_CACHE")
        result = self.run_cache(release="116")
        self.assertNotEqual(result.returncode, 0)
        self.assertIn("changed after failure", old.read_text())
        self.assertTrue((self.directory / "envs/priva_acmg/.priva-vep-release-transition.json").exists())

    def test_direct_helper_cannot_install_cache_for_another_executable(self):
        self.manager("micromamba")
        result = self.run_cache(release="116")
        self.assertNotEqual(result.returncode, 0)
        self.assertIn("Requested VEP 116; executable=113", result.stderr)
        self.assertFalse(any(call[0] == "vep_install" for call in self.calls()))


class PluginReleaseReceiptTests(unittest.TestCase):
    """Small text fixtures exercise publication; no upstream files or downloads."""

    @classmethod
    def setUpClass(cls):
        import importlib.util
        spec = importlib.util.spec_from_file_location("vep_plugins_release", ROOT / "scripts/vep_plugins_release.py")
        cls.helper = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(cls.helper)

    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory(prefix="priva_plugin_release_", dir=ROOT / "tests")
        self.addCleanup(self.temporary.cleanup)
        self.directory = Path(self.temporary.name)
        self.plugins = self.directory / "plugins"
        self.plugins.mkdir()
        self.names = ["SpliceAI", "SpliceVault"]

    def stage(self, release):
        directory = self.directory / ("stage-" + release)
        directory.mkdir(exist_ok=True)
        (directory / "plugin_config.txt").write_text("official config fixture " + release)
        for name in self.names:
            (directory / (name + ".pm")).write_text("package " + name + ";\n# fixture " + release + "\n1;\n")
        return directory

    def test_transition_uses_the_same_lexical_cache_root_normalization(self):
        actual = self.directory / "actual/deep"
        actual.mkdir(parents=True)
        link = self.directory / "link"
        link.symlink_to(actual, target_is_directory=True)
        plugins = self.directory / "cache/CustomPlugins"
        plugins.mkdir(parents=True)
        (plugins / "SpliceAI.pm").write_text("package SpliceAI;\n1;\n")
        prefix = self.directory / "environment"
        prefix.mkdir()
        transition = prefix / ".priva-vep-release-transition.json"
        self.helper.capture_transition(transition, {
            "cache_root": str(link) + "/../cache",
            "cache_folders": {"vep_plugin": "CustomPlugins"},
        }, "113", "116")
        record = json.loads(transition.read_text())
        self.assertEqual(record["plugins_directory"], str(plugins.resolve()))
        self.assertEqual(set(record["plugins"]), {"SpliceAI.pm"})

    def test_identical_unrecorded_files_are_verified_without_replacement(self):
        stage = self.stage("113")
        for name in self.names:
            shutil.copyfile(stage / (name + ".pm"), self.plugins / (name + ".pm"))
        before = {name: (self.plugins / (name + ".pm")).stat().st_ino for name in self.names}
        self.helper.install(self.plugins, "113", self.names, stage)
        self.assertTrue(self.helper.matches(self.plugins, self.helper.load_receipt(self.plugins), "113", self.names))
        self.assertEqual(before, {name: (self.plugins / (name + ".pm")).stat().st_ino for name in self.names})

    def test_html_or_wrong_package_cannot_publish_receipt_or_files(self):
        stage = self.stage("116")
        for text in ("<html>download temporarily unavailable</html>", "package WrongName;\n1;\n"):
            with self.subTest(text=text):
                (stage / "SpliceAI.pm").write_text(text)
                with self.assertRaisesRegex(ValueError, "does not declare package SpliceAI"):
                    self.helper.install(self.plugins, "116", self.names, stage)
                self.assertEqual(list(self.plugins.iterdir()), [])

    def test_partial_release_publication_can_retry_without_claiming_success(self):
        self.helper.install(self.plugins, "113", self.names, self.stage("113"))
        stage = self.stage("116")
        original_replace = self.helper.os.replace

        def fail_second(source, destination):
            if Path(destination).name == "SpliceVault.pm":
                raise OSError("fixture write failure")
            return original_replace(source, destination)

        with mock.patch.object(self.helper.os, "replace", side_effect=fail_second):
            with self.assertRaisesRegex(OSError, "fixture write failure"):
                self.helper.install(self.plugins, "116", self.names, stage)
        receipt = self.helper.load_receipt(self.plugins)
        self.assertEqual(receipt["release"], "113")
        self.assertFalse(self.helper.matches(self.plugins, receipt, "116", self.names))
        self.helper.install(self.plugins, "116", self.names, stage)
        self.assertTrue(self.helper.matches(self.plugins, self.helper.load_receipt(self.plugins), "116", self.names))
        self.assertFalse(list(self.plugins.glob(".priva-plugin-*")))

    def test_partial_retry_still_preserves_unknown_customized_bytes(self):
        self.helper.install(self.plugins, "113", self.names, self.stage("113"))
        stage = self.stage("116")
        shutil.copyfile(stage / "SpliceAI.pm", self.plugins / "SpliceAI.pm")
        (self.plugins / "SpliceVault.pm").write_text("package SpliceVault;\n# customized\n1;\n")
        with self.assertRaisesRegex(ValueError, "differs from official"):
            self.helper.install(self.plugins, "116", self.names, stage)
        self.assertIn("customized", (self.plugins / "SpliceVault.pm").read_text())
        self.assertEqual(self.helper.load_receipt(self.plugins)["release"], "113")


class CachePathTests(unittest.TestCase):
    """Exercise cache configuration entirely in memory; no test folders/files."""

    @classmethod
    def setUpClass(cls):
        source = INSTALLER.read_text()
        program = source.split("<<'PRIVA_CACHE_PY'\n", 1)[1].split("\nPRIVA_CACHE_PY", 1)[0]
        cls.namespace = {"__name__": "priva_cache_tests", "CONFIG_PATHS_HELPER": ROOT / "scripts/installer_config_paths.py"}
        exec(compile(program, str(INSTALLER) + ":cache_paths", "exec"), cls.namespace)

    def resolve(self, config):
        return self.namespace["resolve_cache_paths"](config)

    def test_defaults_need_only_a_root(self):
        directories, updates = self.resolve({"cache_root": "/scratch/caches/"})
        self.assertEqual(directories, {
            "vep_cache_dir": "/scratch/caches/VEP_cache",
            "vep_plugins_dir": "/scratch/caches/VEP_plugin",
            "vep_plugins_cachedir": "/scratch/caches/VEP_plugin_cache",
            "gnomad_vcf_dir": "/scratch/caches/gnomAD",
            "clinvar_vcf_dir": "/scratch/caches/ClinVar",
            "cadd_parent_dir": "/scratch/caches/CADD",
        })
        self.assertEqual(updates, directories)

    def test_partial_overrides_and_null_mapping(self):
        directories, _ = self.resolve({
            "cache_root": "/scratch/caches",
            "cache_folders": {"vep_cache": "VEP_113", "cadd": "CADD_17"},
        })
        self.assertEqual(directories["vep_cache_dir"], "/scratch/caches/VEP_113")
        self.assertEqual(directories["cadd_parent_dir"], "/scratch/caches/CADD_17")
        self.assertEqual(directories["gnomad_vcf_dir"], "/scratch/caches/gnomAD")
        self.assertEqual(self.resolve({"cache_root": "/scratch/caches", "cache_folders": None}),
                         self.resolve({"cache_root": "/scratch/caches"}))

    def test_legacy_configurations_keep_their_full_paths(self):
        self.assertEqual(self.resolve({"vep_cache_dir": "/shared/vep"}), ({}, {}))

    def test_root_change_plans_managed_children_and_external_cache_destinations(self):
        config = {
            "cache_root": "/new/caches", "vep_plugins_cachedir": "/old/plugins",
            "spliceai_snv_prescore": "/old/plugins/SpliceAI/snv.vcf.gz",
            "spliceai_indel_prescore": "/shared/SpliceAI/indel.vcf.gz",
            "primateai_prescore": "/old/plugins_extra/scores.tsv.bgz",
            "input_vcf": "/old/plugins/cohort.vcf.gz",
        }
        _, updates = self.resolve(config)
        self.assertEqual(updates["spliceai_snv_prescore"],
                         "/new/caches/VEP_plugin_cache/SpliceAI/snv.vcf.gz")
        self.assertTrue(updates["spliceai_indel_prescore"].startswith("/new/caches/VEP_plugin_cache/SpliceAI/"))
        self.assertTrue(updates["primateai_prescore"].startswith("/new/caches/VEP_plugin_cache/PrimateAI/"))
        self.assertNotIn("input_vcf", updates)

    def test_folder_rename_and_nested_old_folders(self):
        _, updates = self.resolve({
            "cache_root": "/cache", "cache_folders": {"vep_plugin": "Plugins"},
            "vep_cache_dir": "/cache/VEP", "vep_plugins_dir": "/cache/VEP/plugins",
            "loftee_repo": "/cache/VEP/plugins/LoFtee_repos/loftee-hg38",
        })
        self.assertEqual(updates["loftee_repo"], "/cache/Plugins/LoFtee_repos/loftee-hg38")

    def test_invalid_roots_are_rejected(self):
        for root in ("relative", "/", "/cache root", "/cache;$USER", 12, False, []):
            with self.subTest(root=root), self.assertRaises(ValueError):
                self.resolve({"cache_root": root})

    def test_folder_overrides_require_known_names_and_single_components(self):
        for folders in ([], False, "VEP", {"typo": "VEP"}, {"vep_cache": "../VEP"},
                        {"vep_cache": "/VEP"}, {"vep_cache": "one/two"},
                        {"vep_cache": ""}, {"vep_cache": None}, {"vep_cache": 5}):
            with self.subTest(folders=folders), self.assertRaises(ValueError):
                self.resolve({"cache_root": "/cache", "cache_folders": folders})
        with self.assertRaisesRegex(ValueError, "Indent"):
            self.resolve({"cache_root": "/cache", "vep_cache": "VEP"})

    def test_render_preserves_comments_and_unrelated_values(self):
        original = '# Header\nvep_cache_dir: "/old/#vep"  # keep comment\ninput_vcf: "keep.vcf.gz"\n'
        render = self.namespace["render_cache_config"]
        updated = render(original, {"vep_cache_dir": "/cache/VEP", "cadd_parent_dir": "/cache/CADD"})
        self.assertIn('vep_cache_dir: "/cache/VEP"  # keep comment\n', updated)
        self.assertIn('input_vcf: "keep.vcf.gz"\n', updated)
        self.assertTrue(updated.startswith('# Header\n'))
        self.assertEqual(render(updated, {"vep_cache_dir": "/cache/VEP", "cadd_parent_dir": "/cache/CADD"}), updated)
        with self.assertRaisesRegex(ValueError, "Duplicate"):
            render('vep_cache_dir: a\nvep_cache_dir: b\n', {"vep_cache_dir": "/cache/VEP"})

    def test_check_mode_never_writes_or_creates_directories(self):
        output = io.StringIO()
        with mock.patch.object(Path, "read_text", side_effect=AssertionError("Unexpected read")), \
             mock.patch.object(Path, "write_text", side_effect=AssertionError("Unexpected write")), \
             mock.patch.object(Path, "mkdir", side_effect=AssertionError("Unexpected directory")), \
             contextlib.redirect_stdout(output):
            self.namespace["configure_cache_paths"]("unused", {"cache_root": "/cache"}, "--check")
        self.assertEqual(json.loads(output.getvalue())["gnomad_vcf_dir"], "/cache/gnomAD")

    def test_write_is_idempotent_and_only_creates_selected_directories(self):
        config = {"cache_root": "/cache"}
        directories, updates = self.resolve(config)
        original = self.namespace["render_cache_config"]('cache_root: /cache\n', updates)
        with mock.patch.object(Path, "read_text", return_value=original), \
             mock.patch.object(Path, "write_text") as write, \
             mock.patch.object(Path, "mkdir", autospec=True) as mkdir, \
             contextlib.redirect_stdout(io.StringIO()):
            self.namespace["configure_cache_paths"]("unused", config, "write")
        write.assert_not_called()
        self.assertEqual({str(call.args[0]) for call in mkdir.call_args_list}, set(directories.values()))

    def test_main_stops_if_cache_configuration_fails(self):
        script = '''
            source "$1" >/dev/null || exit
            priva_install_basic_paths() { PRIVA_INSTALL_TMPDIR="$TMPDIR"; }
            conda_install_vep() { PRIVA_INSTALL_ENV_NAME=priva_acmg; }
            # This boundary fixture tests failure in cache-root planning only.
            priva_manual_cache_inputs() { return 0; }
            priva_install_prerequisites() { return 0; }
            nm_to_enst_mapping_install() { return 0; }
            priva_configure_cache_paths() { printf 'CACHE_CHECK_FAILED\\n' >&2; return 1; }
            reference_genome_install() { printf 'UNEXPECTED_REFERENCE\\n'; }
            main_install "$2"
        '''
        result = subprocess.run(["bash", "--noprofile", "--norc", "-c", script,
                                 "bash", str(INSTALLER), str(ROOT / "config.yaml")],
                                capture_output=True, text=True, timeout=10)
        self.assertNotEqual(result.returncode, 0)
        self.assertIn("CACHE_CHECK_FAILED", result.stderr)
        self.assertNotIn("UNEXPECTED_REFERENCE", result.stdout)

    def test_missing_or_empty_root_does_not_require_python(self):
        script = '''
            source "$1" >/dev/null 2>&1 || exit
            yq() { printf '%s\\n' "$PRIVA_TEST_ROOT_VALUE"; }
            python() { printf 'UNEXPECTED_PYTHON\\n'; return 1; }
            priva_configure_cache_paths unused
        '''
        for value in ("null", '\"\"'):
            result = subprocess.run(["bash", "--noprofile", "--norc", "-c", script,
                                     "bash", str(INSTALLER)],
                                    env=dict(os.environ, PRIVA_TEST_ROOT_VALUE=value),
                                    capture_output=True, text=True, timeout=10)
            self.assertEqual(result.returncode, 0, result.stderr)
            self.assertEqual(result.stdout, "")

    def test_yaml_read_error_does_not_fall_back_to_legacy_paths(self):
        script = '''
            source "$1" >/dev/null 2>&1 || exit
            yq() { return 2; }
            priva_configure_cache_paths unused
        '''
        result = subprocess.run(["bash", "--noprofile", "--norc", "-c", script,
                                 "bash", str(INSTALLER)],
                                capture_output=True, text=True, timeout=10)
        self.assertNotEqual(result.returncode, 0)

    def test_cadd_uses_configured_parent_without_prompt_or_environment_search(self):
        script = '''
            source "$1" >/dev/null || exit
            priva_configure_cache_paths() { return 0; }
            read_yaml() {
                case "$2" in
                    cache_root) printf '/scratch/caches' ;;
                    cadd_parent_dir) printf '/scratch/caches/Custom_CADD' ;;
                    cadd_version) printf 'v1.7' ;;
                    cadd_zip_download_url) printf 'https://github.com/kircherlab/CADD-scripts/archive/refs/tags/v1.7.2.zip' ;;
                esac
            }
            read() { printf 'UNEXPECTED_PROMPT\\n' >&2; return 90; }
            find() { printf 'FIND %s\\n' "$1" >&2; }
            python() { printf 'SOURCE_HELPER %s\\n' "$*"; return 73; }
            wget() { printf 'UNEXPECTED_ZIP_DOWNLOAD\\n'; return 90; }
            unzip() { printf 'UNEXPECTED_UNZIP\\n'; return 90; }
            CADD_install unused hg38
        '''
        result = subprocess.run(["bash", "--noprofile", "--norc", "-c", script,
                                 "bash", str(INSTALLER)], capture_output=True, text=True, timeout=10)
        self.assertNotIn("UNEXPECTED_PROMPT", result.stdout + result.stderr)
        searches = [line for line in result.stderr.splitlines() if line.startswith("FIND ")]
        self.assertEqual(searches, ["FIND /scratch/caches/Custom_CADD/"])
        self.assertEqual(result.stdout.strip(),
                         f"SOURCE_HELPER {INSTALLER.parent / 'cadd_runtime_install.py'} "
                         "--config unused --source-only --source "
                         "/scratch/caches/Custom_CADD/CADD-scripts-1.7.2")
        self.assertNotEqual(result.returncode, 0)
        self.assertNotIn("UNEXPECTED_ZIP_DOWNLOAD", result.stdout + result.stderr)
        self.assertNotIn("UNEXPECTED_UNZIP", result.stdout + result.stderr)


if __name__ == "__main__":
    unittest.main()
