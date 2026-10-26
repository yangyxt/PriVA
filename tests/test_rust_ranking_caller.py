"""Exercise optional ranking dispatch and executable identity, without model inference.

The AST pins protect the existing ACMG calculations and write order. Temporary
files and the executable below are dispatch fixtures, not scientific inputs.
"""
import argparse
import ast
import contextlib
import copy
import hashlib
import importlib.util
import json
import os
from pathlib import Path
import shlex
import subprocess
import sys
import tempfile
from types import SimpleNamespace
import unittest
from unittest import mock

ROOT = Path(__file__).resolve().parents[1]
spec = importlib.util.spec_from_file_location("ranking_caller_task", ROOT / "bin/priva_task.py")
adapter = importlib.util.module_from_spec(spec)
spec.loader.exec_module(adapter)


def ast_digest(node):
    return hashlib.sha256(ast.dump(node, include_attributes=False).encode()).hexdigest()


def acmg_tree():
    return ast.parse((ROOT / "scripts/acmg_criteria_assign.py").read_text())


def ranking_dispatch(function):
    candidates = [node for node in function.body if isinstance(node, ast.If)
                  and ast.dump(node.test) == ast.dump(ast.parse("ranking_binary is None", mode="eval").body)]
    if len(candidates) != 1:
        raise AssertionError("Expected one explicit default/native ranking boundary")
    return candidates[0]


class RankingCallBoundary(unittest.TestCase):
    def test_existing_criteria_calculation_and_write_order_are_unchanged(self):
        tree = acmg_tree()
        functions = [node for node in tree.body if isinstance(node, ast.FunctionDef)
                     and node.name == "ACMG_criteria_assign"]
        self.assertEqual(len(functions), 1)
        function = copy.deepcopy(functions[0])
        self.assertEqual(function.args.args[-1].arg, "ranking_binary")
        self.assertEqual(ast.dump(function.args.defaults[-1]), "Constant(value=None)")
        function.args.args.pop()
        function.args.defaults.pop()
        dispatch = ranking_dispatch(function)
        self.assertEqual(len(dispatch.body), 1)
        self.assertEqual(len(dispatch.orelse), 1)
        original_call = dispatch.body[0]
        native_call = copy.deepcopy(dispatch.orelse[0])
        native_keywords = native_call.value.keywords
        self.assertEqual(native_keywords[-1].arg, "ranking_binary")
        self.assertEqual(ast.dump(native_keywords[-1].value), "Name(id='ranking_binary', ctx=Load())")
        native_keywords.pop()
        self.assertEqual(ast.dump(native_call), ast.dump(original_call))
        function.body[function.body.index(dispatch):function.body.index(dispatch) + 1] = dispatch.body
        self.assertEqual(ast_digest(function), "fba4729cf524acc543df0d45a493007d7dc7f1e39ffb1bb8be4c474fac8330be")

    def test_cli_preserves_all_existing_arguments(self):
        main = copy.deepcopy(next(node for node in acmg_tree().body if isinstance(node, ast.If)
                                  and isinstance(node.test, ast.Compare)
                                  and ast.unparse(node.test) == "__name__ == '__main__'"))
        added = [node for node in main.body if isinstance(node, ast.Expr)
                 and isinstance(node.value, ast.Call) and node.value.args
                 and isinstance(node.value.args[0], ast.Constant)
                 and node.value.args[0].value == "--ranking_binary"]
        self.assertEqual(len(added), 1)
        main.body.remove(added[0])
        call = main.body[-1].value
        self.assertEqual(call.func.id, "ACMG_criteria_assign")
        self.assertEqual(call.keywords[-1].arg, "ranking_binary")
        call.keywords.pop()
        self.assertEqual(ast_digest(main), "c1951401d8df01b7fed0c6d7e389563492be36378da909d593386767b42f2e1b")

    def test_actual_argument_parser_forwards_default_and_selected_binary(self):
        main = next(node for node in acmg_tree().body if isinstance(node, ast.If)
                    and isinstance(node.test, ast.Compare)
                    and ast.unparse(node.test) == "__name__ == '__main__'")
        required = []
        for node in main.body:
            if isinstance(node, ast.Expr) and isinstance(node.value, ast.Call):
                if any(keyword.arg == "required" and isinstance(keyword.value, ast.Constant)
                       and keyword.value.value is True for keyword in node.value.keywords):
                    required.extend([node.value.args[0].value, "fixture"])
        observed = []
        for binary in (None, "/qualified/bin/priva"):
            calls = []
            def capture(*args, **kwargs):
                calls.append((args, kwargs))
                return None, None
            namespace = {"__name__": "__main__", "ap": argparse,
                         "ACMG_criteria_assign": capture,
                         "DEFAULT_DISPENSABLE_GENE_LIST": "bundled-genes.txt",
                         "DEFAULT_HPO_CONDITION_MECHANISM_CACHE": Path("conditions.json"),
                         "DEFAULT_MECHANISM_JSON": Path("mechanisms.json")}
            argv = ["acmg_criteria_assign.py", *required]
            if binary is not None:
                argv.extend(["--ranking_binary", binary])
            with mock.patch.object(sys, "argv", argv):
                exec(compile(ast.Module(body=[copy.deepcopy(main)], type_ignores=[]), "actual-cli", "exec"), namespace)
            self.assertEqual(len(calls), 1)
            self.assertEqual(calls[0][1].pop("ranking_binary"), binary)
            observed.append(calls[0])
        self.assertEqual(observed[0], observed[1])

    def test_actual_rank_call_uses_current_frame_once(self):
        function = next(node for node in acmg_tree().body if isinstance(node, ast.FunctionDef)
                        and node.name == "ACMG_criteria_assign")
        boundary = ranking_dispatch(function)
        for binary in (None, "/qualified/bin/priva"):
            calls = []
            result = object()
            def capture(*args, **kwargs):
                calls.append((args, kwargs))
                return result
            variables = {node.id: object() for node in ast.walk(boundary) if isinstance(node, ast.Name)}
            original = dict(variables)
            variables.update(sort_and_rank_variants=capture, ranking_binary=binary)
            exec(compile(ast.Module(body=[copy.deepcopy(boundary)], type_ignores=[]), "actual-boundary", "exec"), variables)
            self.assertEqual(len(calls), 1)
            self.assertIs(variables["anno_df"], result)
            self.assertEqual(calls[0][0], tuple(original[name] for name in ("anno_df", "ped_df", "fam_name")))
            for key, value in calls[0][1].items():
                self.assertIs(value, binary if key == "ranking_binary" else original[key])
            self.assertEqual("ranking_binary" in calls[0][1], binary is not None)


class RankingExecutableIdentity(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.folder = Path(self.temp.name)
        self.repo = self.folder / "source"
        self.repo.mkdir()
        for name in ("scripts", "bin", "modules", "conf"):
            (self.repo / name).mkdir()
        for name in ("main.nf", "nextflow.config"):
            (self.repo / name).write_text("fixture\n")
        self.config = {"assembly": "hg19", "threads": 1, "threads_per_fam": 1,
                       "clinvar_gene_stat": "gene.pkl", "repeat_region_file_name": "repeats.bed"}
        for key in adapter.PRIORITY_KEYS:
            path = self.folder / (key + ".fixture")
            path.write_text("dispatch fixture\n")
            self.config[key] = str(path)
        for name in ("data/hpo/genes_to_phenotype.assertions.tsv.gz",
                     "data/alphamissense/alphamissense_mean_score.tsv", "data/ClinVar/gene.pkl",
                     "data/repeats/repeats.bed", "data/hgnc/non_alt_loci_set.tsv",
                     "data/dispensable_genes/dispensable_gene_list.txt"):
            path = self.repo / name
            path.parent.mkdir(parents=True, exist_ok=True)
            path.write_text("dispatch fixture\n")
        self.binary = self.folder / "fixture-ranker"
        self.binary.write_bytes(b"#!/bin/sh\nexit 0\n")
        self.binary.chmod(0o555)
        self.config_path = self.folder / "config.yaml"
        self.manifest_path = self.folder / "resources.json"

    def manifest(self, selected=True, checksum=False):
        if selected:
            self.config["priva_ranking_binary"] = str(self.binary)
        else:
            self.config.pop("priva_ranking_binary", None)
        self.config_path.write_text(adapter.yaml.safe_dump(self.config))
        adapter.build_manifest(SimpleNamespace(config=str(self.config_path), legacy_root=str(self.repo),
            mode="annotated", vcf_hub_mode="off", cadd_hub_mode="copy", checksum=checksum,
            output=str(self.manifest_path)))
        return json.loads(self.manifest_path.read_text())

    def verify(self, manifest):
        return adapter.verify_manifest(manifest, self.config_path, self.repo, "annotated", "off", "copy")

    def alter_bytes_preserving_size_mtime_mode(self):
        before = self.binary.stat()
        self.binary.chmod(0o755)
        self.binary.write_bytes(b"#!/bin/sh\nexit 1\n")
        self.binary.chmod(0o555)
        os.utime(self.binary, ns=(before.st_atime_ns, before.st_mtime_ns))
        self.assertEqual(self.binary.stat().st_size, before.st_size)
        self.assertEqual(self.binary.stat().st_mtime_ns, before.st_mtime_ns)
        self.assertEqual(self.binary.stat().st_mode, before.st_mode)

    def test_default_inventory_and_optional_content_pin(self):
        default = self.manifest(selected=False)
        self.assertEqual(default["identity_kind"], "size-mtime-inventory")
        self.assertNotIn("priva_ranking_binary", self.verify(default))
        selected = self.manifest()
        self.assertEqual(selected["identity_kind"], "size-mtime-inventory+ranking-sha256")
        binary = selected["resources"]["priva_ranking_binary"]
        self.assertEqual(binary["sha256"], hashlib.sha256(self.binary.read_bytes()).hexdigest())
        self.assertEqual(binary["executable_mode"], 0o555)
        self.assertNotIn("sha256", selected["resources"][adapter.PRIORITY_KEYS[0]])
        self.verify(selected)
        self.alter_bytes_preserving_size_mtime_mode()
        with self.assertRaisesRegex(adapter.ContractError, "Resource changed.*priva_ranking_binary"):
            self.verify(selected)

    def test_invalid_binary_is_not_replaced_or_ignored(self):
        paths = (self.folder / "missing", self.folder / "empty", self.folder / "not-executable", self.folder / "directory")
        paths[1].touch()
        paths[1].chmod(0o555)
        paths[2].write_text("nonexecutable fixture\n")
        paths[2].chmod(0o444)
        paths[3].mkdir()
        for path in paths:
            with self.subTest(path=path.name):
                self.config["priva_ranking_binary"] = str(path)
                with self.assertRaises(adapter.ContractError):
                    adapter.runtime_resources(self.config, self.repo, "annotated")

    def test_changed_during_checksum_capture_fails(self):
        original = adapter.sha256
        def change(path):
            digest = original(path)
            if Path(path) == self.binary:
                self.alter_bytes_preserving_size_mtime_mode()
            return digest
        with mock.patch.object(adapter, "sha256", side_effect=change):
            with self.assertRaisesRegex(adapter.ContractError, "changed during identity capture"):
                adapter.ranking_binary_record(self.binary)

    def test_relaunch_checks_binary_before_cached_run_identity(self):
        self.manifest()
        inputs = []
        for name in ("input.vcf.gz", "input.vcf.gz.tbi", "family.ped", "paired-cadd.tsv"):
            path = self.folder / name
            path.write_text("launch identity fixture; no scientific execution\n")
            inputs.append(str(path))
        args = SimpleNamespace(mode="annotated", stub=False, legacy_root=str(self.repo),
            source_digest=adapter.source_digest(self.repo), config=str(self.config_path), assembly="hg19",
            annotation_cpus=1, family_cpus=1, resources=str(self.manifest_path), software_id="fixture",
            family="FAM", ped=inputs[2], vcf=inputs[0], vcf_index=inputs[1], cadd=inputs[3],
            vcf_hub_mode="off", cadd_hub_mode="copy", task_boundary="off", outdir=str(self.folder / "output"))
        adapter.launch(args)
        marker = Path(args.outdir) / ".priva-run.json"
        before = marker.read_bytes()
        adapter.launch(args)
        self.assertEqual(marker.read_bytes(), before)
        self.alter_bytes_preserving_size_mtime_mode()
        with self.assertRaisesRegex(adapter.ContractError, "Resource changed.*priva_ranking_binary"):
            adapter.launch(args)
        self.assertEqual(marker.read_bytes(), before)

    def test_private_task_keeps_selected_executable_shared_and_readonly(self):
        self.manifest()
        bundle = self.folder / "bundle"
        bundle.mkdir()
        (bundle / "legacy.yaml").write_bytes(self.config_path.read_bytes())
        (bundle / "family.ped").write_text("FAM\tchild\t0\t0\t0\t2\n")
        task = self.folder / "task"
        task.mkdir()
        before = self.binary.read_bytes(), self.binary.stat().st_mode
        with contextlib.chdir(task):
            config_path, _ = adapter.private_config(bundle, {"base": "input", "assembly": "hg19", "mode": "annotated"}, "prioritization", 1)
            actual = adapter.read_config(config_path)
        self.assertEqual(actual["priva_ranking_binary"], str(self.binary))
        self.assertFalse((task / "result/private_inputs/priva_ranking_binary").exists())
        self.assertEqual((self.binary.read_bytes(), self.binary.stat().st_mode), before)

    def test_scientific_qualification_gates_remain_closed(self):
        tree = ast.parse((ROOT / "bin/priva_task.py").read_text())
        expected = {"require_qualified_execution": "6d4bd61d7fbb54b107005b4d97978257bb83de2a93f0238718d19b12e2ce6cdd",
                    "read_config": "e675eff9923a92114ae12f1b13b045bcd20f54bb9e78dae8498e0df58877959b"}
        for node in tree.body:
            if isinstance(node, ast.FunctionDef) and node.name in expected:
                self.assertEqual(ast_digest(node), expected.pop(node.name))
        self.assertFalse(expected)
        with self.assertRaisesRegex(adapter.ContractError, "Real raw annotation is gated"):
            adapter.require_qualified_execution("raw", False)


class RankingShellDispatch(unittest.TestCase):
    def test_actual_shell_preserves_default_skip_and_passes_explicit_binary_once(self):
        source = (ROOT / "scripts/prioritization_vcf_per_fam.sh").read_text()
        function = source.split("function assign_acmg_criteria () {", 1)[1].split("\nfunction ", 1)[0]
        function = "function assign_acmg_criteria () {" + function
        with tempfile.TemporaryDirectory() as folder:
            root = Path(folder)
            old = root / "old.fixture"
            old.write_text("fixture\n")
            for name in ("data/alphamissense/alphamissense_mean_score.tsv", "data/repeats/repeats.bed",
                         "scripts/acmg_criteria_assign.py", "scripts/acmg_pm1_pp2_bp1_missense.py"):
                path = root / name
                path.parent.mkdir(parents=True, exist_ok=True)
                path.write_text("fixture\n")
                os.utime(path, ns=(1_000_000_000, 1_000_000_000))
            os.utime(old, ns=(1_000_000_000, 1_000_000_000))
            input_table = root / "input.tsv"
            input_table.write_text("fixture\n")
            os.utime(input_table, ns=(2_000_000_000, 2_000_000_000))
            output = root / "input.acmg.tsv"
            output.write_text("fixture\n")
            os.utime(output, ns=(3_000_000_000, 3_000_000_000))
            binary = root / "fixture-ranker"
            binary.write_text("#!/bin/sh\nexit 0\n")
            binary.chmod(0o555)
            for selected, expected_return, expected_call in (("null", 0, False), (str(binary), 0, True), (str(root / "missing"), 1, False)):
                calls = root / "argv.bin"
                events = root / "events.txt"
                calls.unlink(missing_ok=True)
                events.unlink(missing_ok=True)
                prelude = "\n".join([
                    "BASE_DIR=" + shlex.quote(str(root)), "DATA_DIR=$BASE_DIR/data", "SCRIPT_DIR=$BASE_DIR/scripts",
                    "OLD=" + shlex.quote(str(old)), "CALLS=" + shlex.quote(str(calls)), "EVENTS=" + shlex.quote(str(events)),
                    "SELECTED=" + shlex.quote(selected),
                    'read_yaml() { case "$2" in priva_ranking_binary) printf "%s\\n" "$SELECTED";; assembly) echo hg19;; repeat_region_file_name) echo repeats.bed;; pext_tissues|pp1_vcf|pp1_ped|alt_disease_vcf) echo null;; *) printf "%s\\n" "$OLD";; esac; }',
                    'interpret_splicing_annotations() { echo splicing >> "$EVENTS"; }',
                    'interpret_utr_annotations() { echo utr >> "$EVENTS"; }',
                    'check_path() { return 0; }', 'check_table_column() { return 0; }',
                    'check_vcf_validity() { return 1; }', 'log() { :; }', 'display_table() { :; }',
                    'python() { printf "%s\\0" "$@" > "$CALLS"; echo python >> "$EVENTS"; }', function,
                    "assign_acmg_criteria -i " + shlex.quote(str(input_table)) + " -c " + shlex.quote(str(old)) + " -f FAM -t 1",
                ])
                result = subprocess.run(["bash", "-s"], input=prelude, capture_output=True, text=True, timeout=30)
                self.assertEqual(result.returncode, expected_return, result.stderr)
                self.assertEqual(calls.exists(), expected_call)
                self.assertEqual(events.read_text().splitlines(), ["splicing", "utr"] + (["python"] if expected_call else []))
                if expected_call:
                    argv = calls.read_bytes().rstrip(b"\0").decode().split("\0")
                    self.assertEqual(argv.count("--ranking_binary"), 1)
                    self.assertEqual(argv[argv.index("--ranking_binary") + 1], selected)
                    self.assertEqual(argv[0], str(root / "scripts/acmg_criteria_assign.py"))


if __name__ == "__main__":
    unittest.main()
