"""Exercise the runtime's actual gene-stat path selection; scientific assignment is stubbed."""
import gzip
import json
import os
from pathlib import Path
import pickle
import re
import subprocess
import sys
import tempfile
import unittest


ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT / "scripts/prioritization_vcf_per_fam.sh"
PRELUDE = r'''
set -eo pipefail
log() { printf '%s\n' "$*" >&2; }
read_yaml() {
    case "$2" in
        clinvar_gene_stat) printf '%s\n' "$TEST_STAT" ;;
        assembly) printf 'hg19\n' ;;
        repeat_region_file_name) printf 'resource.txt\n' ;;
        apply_clingen_override) printf 'true\n' ;;
        extreme_rare_PAF|exp_disease_incidence) printf '0.0001\n' ;;
        pp1_vcf|pp1_ped|alt_disease_vcf|ped_file|relevant_gene_list|dispensable_gene_list|pext_tissues) printf '\n' ;;
        *) printf '%s\n' "$TEST_RESOURCE" ;;
    esac
}
interpret_splicing_annotations() { :; }
interpret_utr_annotations() { :; }
check_path() { [[ $# -eq 3 && "$2" == file && -s "$1" ]]; }
check_vcf_validity() { return 1; }
check_table_column() { return 1; }
display_table() { :; }
python() {
    "$TEST_PYTHON" - "$@" <<'ASSIGNMENT_STUB'
import gzip, json, os, pickle, sys
arguments = sys.argv[1:]
selected = arguments[arguments.index("--clinvar_gene_stat_pkl") + 1]
with gzip.open(selected, "rb") as handle:
    counts = pickle.load(handle)
with open(os.environ["TEST_EVENTS"], "a") as handle:
    handle.write(json.dumps({"path": selected, "counts": counts, "arguments": arguments}) + "\n")
ASSIGNMENT_STUB
}
'''


class ClinVarRuntimePaths(unittest.TestCase):
    def setUp(self):
        temporary = tempfile.TemporaryDirectory(
            prefix="priva_clinvar_runtime_", dir=os.environ.get("PRIVA_TEST_TMPDIR", os.environ.get("TMPDIR")))
        self.addCleanup(temporary.cleanup)
        self.folder = Path(temporary.name)
        self.data = self.folder / "source/data"
        self.baseline = self.data / "ClinVar/clinvar_2star_stats.pkl.gz"
        self.write_stats(self.baseline, {"GENE": {"clinvar_pathogenic": 1}})
        self.selected = self.folder / "fresh statistics's folder/clinvar_2star_stats.pkl.gz"
        self.fresh_counts = {"GENE": {"clinvar_pathogenic": 17}, "OTHER": {"clinvar_pathogenic": 3}}
        self.write_stats(self.selected, self.fresh_counts)
        self.resource = self.folder / "resource.txt"
        self.resource.write_text("unrelated routing fixture\n")
        for relative in ("repeats/resource.txt", "alphamissense/alphamissense_mean_score.tsv"):
            path = self.data / relative
            path.parent.mkdir(parents=True, exist_ok=True)
            path.write_text("unrelated bundled fixture\n")
        self.input = self.folder / "input.tsv"
        self.input.write_text("routing fixture\n")
        self.events = self.folder / "events.jsonl"
        self.events.touch()
        source = SOURCE.read_text()
        self.function = re.search(r"^function assign_acmg_criteria\s*\([^\n]*\)\s*\{.*?^\}",
                                  source, re.M | re.S).group()
        self.env = dict(os.environ, TEST_PYTHON=sys.executable, TEST_EVENTS=str(self.events),
                        TEST_RESOURCE=str(self.resource), TEST_INPUT=str(self.input),
                        BASE_DIR=str(self.folder / "source"), DATA_DIR=str(self.data), SCRIPT_DIR=str(self.folder))

    @staticmethod
    def write_stats(path, counts):
        path.parent.mkdir(parents=True, exist_ok=True)
        with gzip.open(path, "wb") as handle:
            pickle.dump(counts, handle)

    def run_assignment(self, selected):
        env = dict(self.env, TEST_STAT=str(selected))
        result = subprocess.run(["bash", "--noprofile", "--norc", "-c",
                                 PRELUDE + self.function + '\nassign_acmg_criteria -i "$TEST_INPUT" -c fixture.yaml -t 1'],
                                cwd=self.folder, env=env, text=True, capture_output=True, timeout=30)
        events = [json.loads(line) for line in self.events.read_text().splitlines()]
        return result, events

    def test_absolute_selected_path_with_spaces_and_quote_reaches_reader_once(self):
        before = self.baseline.read_bytes()
        result, events = self.run_assignment(self.selected)
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertEqual(len(events), 1)
        self.assertEqual(events[0]["path"], str(self.selected))
        self.assertEqual(events[0]["counts"], self.fresh_counts)
        self.assertEqual(events[0]["arguments"].count("--clinvar_gene_stat_pkl"), 1)
        self.assertEqual(self.baseline.read_bytes(), before)

    def test_relative_path_with_directory_is_preserved_from_caller_directory(self):
        selected = self.selected.relative_to(self.folder)
        result, events = self.run_assignment(selected)
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertEqual(events[0]["path"], str(selected))
        self.assertEqual(events[0]["counts"], self.fresh_counts)

    def test_bare_name_keeps_legacy_bundled_location(self):
        result, events = self.run_assignment(self.baseline.name)
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertEqual(events[0]["path"], str(self.baseline))
        self.assertEqual(events[0]["counts"]["GENE"]["clinvar_pathogenic"], 1)

    def test_bare_name_with_spaces_is_one_bundled_filename(self):
        named = self.baseline.with_name("custom gene statistics.pkl.gz")
        self.write_stats(named, self.fresh_counts)
        result, events = self.run_assignment(named.name)
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertEqual(events[0]["path"], str(named))
        self.assertEqual(events[0]["counts"], self.fresh_counts)

    def test_explicit_missing_selection_never_falls_back_to_existing_same_named_baseline(self):
        self.selected.unlink()
        result, events = self.run_assignment(self.selected)
        self.assertNotEqual(result.returncode, 0)
        self.assertEqual(events, [])
        self.assertTrue(self.baseline.is_file())


if __name__ == "__main__":
    unittest.main()
