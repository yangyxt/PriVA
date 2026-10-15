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
import shlex
import sys

root = Path(os.environ["PRIVA_TEST_ROOT"])
args = sys.argv[1:]
kind = Path(sys.argv[0]).name
with (root / "calls.jsonl").open("a") as handle:
    handle.write(json.dumps([kind, *args]) + "\n")

def option(flag):
    return args[args.index(flag) + 1]

def create_environment(name):
    directory = root / "envs" / name / "bin"
    directory.mkdir(parents=True, exist_ok=True)
    for tool in ("yq", "samtools", "vep"):
        if tool == os.environ.get("PRIVA_TEST_MISSING_TOOL"):
            continue
        destination = directory / tool
        destination.write_text(Path(os.environ["PRIVA_TEST_TOOL_TEMPLATE"]).read_text())
        destination.chmod(0o755)

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
    assert Path(option("--file")).is_file()
    create_environment(option("--name"))
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
    assert (root / "envs" / option("--name")).is_dir()
    for name in ("ensembl-vep", "perl-bio-procedural", "perl-bioperl"):
        print(name, "1.0 test")
elif args[:4] == ["env", "config", "vars", "set"]:
    assert option("--prefix") == os.environ["CONDA_PREFIX"]
    if os.environ.get("PRIVA_TEST_FAIL_VARS"):
        sys.exit(44)
else:
    raise SystemExit("Unexpected fake manager call: " + repr(args))
'''

FAKE_TOOL = r'''
import json
import os
from pathlib import Path
import sys

root = Path(os.environ["PRIVA_TEST_ROOT"])
tool = Path(sys.argv[0]).name
with (root / "calls.jsonl").open("a") as handle:
    handle.write(json.dumps([tool, *sys.argv[1:]]) + "\n")
if tool == "yq":
    assert Path(os.environ["CONDA_PREFIX"]) / "bin" == Path(sys.argv[0]).parent
    config = json.loads(Path(sys.argv[2]).read_text())
    print(json.dumps(config.get(sys.argv[1].lstrip("."))))
'''


class InstallerBootstrapTests(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory(prefix="priva_bootstrap_", dir=ROOT / "tests")
        self.addCleanup(self.temporary.cleanup)
        self.directory = Path(self.temporary.name)
        self.bin = self.directory / "bin"
        self.bin.mkdir()
        (self.directory / "envs").mkdir()
        self.bash = shutil.which("bash")
        for command in ("bash", "sed", "awk", "dirname", "realpath", "date", "basename", "grep", "head"):
            executable = shutil.which(command)
            self.assertIsNotNone(executable, command)
            (self.bin / command).symlink_to(executable)
        self.tool_template = self.directory / "tool_template.py"
        self.write_executable(self.tool_template, FAKE_TOOL)
        self.config = self.directory / "config.json"
        self.config.write_text(json.dumps({
            "assembly": "hg38", "ref_genome": str(self.directory / "reference.fa"),
            "vep_cache_dir": str(self.directory), "vep_plugins_dir": str(self.directory),
            "vep_plugins_cachedir": str(self.directory),
            "conda_env_yaml": "/old/machine/ignored_environment.yml",
        }))
        self.env = {key: value for key, value in os.environ.items()
                    if not key.startswith(("CONDA", "MAMBA", "PRIVA_", "BASH_FUNC_"))
                    and key not in ("BASH_ENV", "ENV")}
        self.env.update({
            "PATH": str(self.bin), "PRIVA_TEST_ROOT": str(self.directory),
            "PRIVA_TEST_TOOL_TEMPLATE": str(self.tool_template),
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
        for tool in ("yq", "samtools", "vep"):
            if tool != missing_tool:
                self.write_executable(directory / tool, FAKE_TOOL)

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
                self.assertEqual(create[create.index("--file") + 1], str(ROOT / "acmg_conda_rough.yml"))
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

    def test_vep_failure_is_not_masked_by_environment_setup(self):
        self.manager("micromamba")
        result = self.run_shell('''
            conda_install_vep || exit
            getopt() { printf -- '--VEP_VERSION 113 --\n'; }
            perl() { printf '%s/lib/site_perl\n' "$CONDA_PREFIX"; }
            vep_install() { return 23; }
            vep_cache_api_install --VEP_VERSION 113
        ''')
        self.assertNotEqual(result.returncode, 0)
        self.assertIn("Failed to install the VEP API", result.stderr)
        self.assertFalse(any(call[1:5] == ["env", "config", "vars", "set"] for call in self.calls()))


class CachePathTests(unittest.TestCase):
    """Exercise cache configuration entirely in memory; no test folders/files."""

    @classmethod
    def setUpClass(cls):
        source = INSTALLER.read_text()
        program = source.split("<<'PRIVA_CACHE_PY'\n", 1)[1].split("\nPRIVA_CACHE_PY", 1)[0]
        cls.namespace = {"__name__": "priva_cache_tests"}
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

    def test_root_change_relocates_children_and_preserves_external_files(self):
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
        for key in ("spliceai_indel_prescore", "primateai_prescore", "input_vcf"):
            self.assertNotIn(key, updates)

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
            conda_install_vep() { PRIVA_INSTALL_ENV_NAME=priva_acmg; }
            priva_configure_cache_paths() { printf 'CACHE_CHECK_FAILED\\n'; return 1; }
            reference_genome_install() { printf 'UNEXPECTED_REFERENCE\\n'; }
            main_install "$2"
        '''
        result = subprocess.run(["bash", "--noprofile", "--norc", "-c", script,
                                 "bash", str(INSTALLER), str(ROOT / "config.yaml")],
                                capture_output=True, text=True, timeout=10)
        self.assertNotEqual(result.returncode, 0)
        self.assertIn("CACHE_CHECK_FAILED", result.stdout)
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
                    cadd_zip_download_url) printf 'https://example.invalid/cadd.zip' ;;
                esac
            }
            read() { printf 'UNEXPECTED_PROMPT\\n' >&2; return 90; }
            find() { printf 'FIND %s\\n' "$1" >&2; }
            wget() { printf 'DOWNLOAD %s\\n' "$*"; return 1; }
            CADD_install unused hg38
        '''
        result = subprocess.run(["bash", "--noprofile", "--norc", "-c", script,
                                 "bash", str(INSTALLER)], capture_output=True, text=True, timeout=10)
        self.assertNotIn("UNEXPECTED_PROMPT", result.stdout + result.stderr)
        searches = [line for line in result.stderr.splitlines() if line.startswith("FIND ")]
        self.assertEqual(searches, ["FIND /scratch/caches/Custom_CADD/"])
        self.assertIn("-O /scratch/caches/Custom_CADD/CADD-scripts.zip", result.stdout)


if __name__ == "__main__":
    unittest.main()
