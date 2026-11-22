"""HPO shell routing and freshness using actual files; conversion is explicitly stubbed."""
import gzip
import json
import os
from pathlib import Path
import re
import subprocess
import sys
import tempfile
import time
import unittest


ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT / "scripts/prioritization_vcf_per_fam.sh"
PRELUDE = r'''
set -eo pipefail
log() { printf '%s\n' "$*" >&2; }
read_yaml() {
    "$TEST_PYTHON" -c 'import json,sys; v=json.load(open(sys.argv[1])).get(sys.argv[2]); print("null" if v is None else v)' "$1" "$2"
}
check_vcf_validity() { [[ -s "$1" ]]; }
check_table_column() { [[ -s "$1" && "$2" == vep_consq_lof ]]; }
display_table() { [[ -s "$1" ]]; }
count_vcf_records() { printf '1\n'; }
python() {
    "$TEST_PYTHON" - "$@" <<'CONVERSION_STUB'
import json, os, pathlib, sys
arguments = sys.argv[1:]
with open(os.environ["TEST_EVENTS"], "a") as stream:
    stream.write(json.dumps(["conversion", *arguments]) + "\n")
pathlib.Path(arguments[arguments.index("-o") + 1]).write_text("vep_consq_lof\nfixture\n")
CONVERSION_STUB
}
assign_acmg_criteria() {
    "$TEST_PYTHON" -c 'import json,os,sys; open(os.environ["TEST_EVENTS"],"a").write(json.dumps(["assignment",*sys.argv[1:]])+"\n")' "$@"
}
'''


class HPORuntimePaths(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory(
            prefix="priva_hpo_paths_", dir=os.environ.get("PRIVA_TEST_TMPDIR", os.environ.get("TMPDIR")))
        self.addCleanup(self.temporary.cleanup)
        self.folder = Path(self.temporary.name)
        self.source = self.folder / "source"
        self.scripts = self.folder / "scripts with ' quote"
        self.scripts.mkdir()
        self.converter = self.scripts / "combine_annotations.py"
        self.converter.write_text("# conversion stub is separate from routing under test\n")
        self.default = self.source / "data/hpo/genes_to_phenotype.assertions.tsv.gz"
        self.default.parent.mkdir(parents=True)
        self.hpo(self.default)
        self.selected = self.folder / "data side's HPO.tsv.gz"
        self.hpo(self.selected)
        self.vcf = self.folder / "public.vcf.gz"
        self.vcf.write_text("public routing fixture\n")
        self.cadd = self.folder / "public.cadd.tsv"
        self.cadd.write_text("#header\npublic routing fixture\n")
        self.output = self.folder / "public.tsv"
        self.config_file = self.folder / "config's path with spaces.json"
        self.settings = {"hpo_assertions": str(self.selected), "threads_per_fam": 1}
        self.save_config()
        self.events = self.folder / "events.jsonl"
        self.events.touch()
        source = SOURCE.read_text()
        self.functions = "\n".join(re.search(
            r"^function " + name + r"\s*\([^\n]*\)\s*\{.*?^\}", source, re.S | re.M).group()
            for name in ("prepare_combined_tab", "main_prioritization"))
        self.env = dict(os.environ, TEST_PYTHON=sys.executable, TEST_EVENTS=str(self.events),
                        TEST_CONFIG=str(self.config_file), TEST_VCF=str(self.vcf), TEST_CADD=str(self.cadd),
                        BASE_DIR=str(self.source), SCRIPT_DIR=str(self.scripts))

    @staticmethod
    def hpo(path):
        with gzip.open(path, "wt") as stream:
            stream.write("gene_symbol\tdisease_id\thpo_id\tfrequency\tevidence\treference\n"
                         "GENE\tOMIM:1\tHP:0000006\t-\tTAS\tPMID:1\n")

    def save_config(self):
        self.config_file.write_text(json.dumps(self.settings))

    def run_shell(self, command=None):
        command = command or 'prepare_combined_tab "$TEST_VCF" "$TEST_CADD" 1 "$TEST_CONFIG"'
        result = subprocess.run(["bash", "--noprofile", "--norc", "-c", PRELUDE + self.functions + "\n" + command],
                                env=self.env, cwd=self.folder, text=True, capture_output=True, timeout=30)
        return result, [json.loads(line) for line in self.events.read_text().splitlines()]

    def current_output(self):
        old = time.time() - 100
        for path in (self.vcf, self.selected, self.default, self.config_file, self.converter):
            os.utime(path, (old, old))
        self.output.write_text("vep_consq_lof\ncached fixture\n")
        os.utime(self.output, (old + 50, old + 50))

    def test_explicit_hpo_with_spaces_and_quotes_is_passed_verbatim_in_original_option_order(self):
        result, events = self.run_shell()
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertEqual(events, [["conversion", str(self.converter), "-i", str(self.vcf), "-c", str(self.cadd),
                                   "-o", str(self.output), "-p", str(self.selected), "-t", "1"]])

    def test_three_argument_call_preserves_legacy_default(self):
        result, events = self.run_shell('prepare_combined_tab "$TEST_VCF" "$TEST_CADD" 1')
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertEqual(events[0][events[0].index("-p") + 1], str(self.default))

    def test_unset_empty_and_null_keys_keep_legacy_default(self):
        for value in (None, "", "null"):
            with self.subTest(value=value):
                self.events.write_text("")
                self.settings["hpo_assertions"] = value
                self.save_config()
                result, events = self.run_shell()
                self.assertEqual(result.returncode, 0, result.stderr)
                self.assertEqual(events[0][events[0].index("-p") + 1], str(self.default))

    def test_explicit_missing_or_empty_hpo_fails_without_using_existing_default_or_cached_table(self):
        self.current_output()
        before = self.output.read_bytes()
        for present in (False, True):
            with self.subTest(present=present):
                if present:
                    self.selected.write_bytes(b"")
                else:
                    self.selected.unlink()
                result, events = self.run_shell()
                self.assertNotEqual(result.returncode, 0)
                self.assertEqual(events, [])
                self.assertIn(str(self.selected), result.stderr)
                self.assertEqual(self.output.read_bytes(), before)
                self.assertTrue(self.default.is_file())

    def test_missing_configuration_is_not_replaced_by_default_hpo(self):
        self.config_file.unlink()
        result, events = self.run_shell()
        self.assertNotEqual(result.returncode, 0)
        self.assertEqual(events, [])

    def test_current_table_reuses_selected_hpo_and_configuration(self):
        self.current_output()
        before = self.output.read_bytes()
        result, events = self.run_shell()
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertEqual(events, [])
        self.assertEqual(self.output.read_bytes(), before)

    def test_selected_hpo_refresh_invalidates_cached_table(self):
        self.current_output()
        os.utime(self.selected, (time.time(), time.time()))
        result, events = self.run_shell()
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertEqual(events[0][events[0].index("-p") + 1], str(self.selected))

    def test_changed_configuration_selecting_older_hpo_invalidates_cached_table(self):
        self.current_output()
        self.settings["hpo_assertions"] = str(self.default)
        self.save_config()
        result, events = self.run_shell()
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertEqual(events[0][events[0].index("-p") + 1], str(self.default))

    def test_unselected_legacy_hpo_refresh_does_not_invalidate_explicit_selection(self):
        self.current_output()
        os.utime(self.default, (time.time(), time.time()))
        result, events = self.run_shell()
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertEqual(events, [])

    def test_main_passes_quoted_configuration_to_combination_before_assignment(self):
        result, events = self.run_shell(
            'main_prioritization -i "$TEST_VCF" -c "$TEST_CONFIG" -f PublicFixture -t "$TEST_CADD"')
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertEqual([event[0] for event in events], ["conversion", "assignment"])
        self.assertEqual(events[0][events[0].index("-p") + 1], str(self.selected))
        self.assertEqual(events[1][events[1].index("-c") + 1], str(self.config_file))


if __name__ == "__main__":
    unittest.main()
