"""VEP worker command/lifecycle fixtures, not real VEP scientific equivalence.

The oracle is the exact frozen v7 serial function. A tiny VEP executable records
argv and copies records; real bcftools/tabix exercise the retained serial tail.
Real official VEP 22-versus-8-fork output comparison remains an adoption gate.
"""
from __future__ import annotations

import gzip
import hashlib
import json
import os
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest

import yaml

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
import clinvar_gene_refresh as refresh
from test_clinvar_cache_candidate import write_vcf

ORACLE = ROOT / "tests/fixtures/vep_serial_v7_5be1b946.sh"
ORACLE_SHA = "78cc9a45f457a0cbcc90be5b3379b0a540b92855a30bef4a398e5f0cff70ad80"
FUNCTION_SHA = "f887aa2298360af1779eb9ad28a0fbd737093c66d73e10658af4950bd385415e"

SHELL = r'''
load_script() { source "$1"; }
load_script "$WORKER_SCRIPT" >/dev/null || exit
[[ "$WORKER_ORACLE" != 1 ]] || source "$WORKER_ORACLE_PATH"
check_path() { return 0; }
check_vcf_infotags() { [[ "${WORKER_EXISTING_CSQ:-0}" == 1 ]]; }
randomID() { printf 'FIXED\n'; }
event() { printf '%s\n' "$1" >> "$WORKER_EVENTS"; }
bcftools() { event sort; command bcftools "$@"; }
tabix() { event index; command tabix "$@"; }
announce_remove_tmps() { event remove; rm -f -- "$@"; }
display_vcf() { event display; }
vep() {
    event vep
    python - "$@" <<'PY_VEP'
import gzip, json, os, pathlib, sys
args = sys.argv[1:]
pathlib.Path(os.environ['WORKER_TRACE']).write_text(json.dumps(args))
source = pathlib.Path(args[args.index('-i') + 1])
payload = gzip.open(source, 'rb').read() if source.suffix == '.gz' else source.read_bytes()
output = pathlib.Path(args[args.index('-o') + 1])
status = int(os.environ.get('WORKER_VEP_STATUS', '0'))
output.write_bytes(payload[:100] if status else payload)
pathlib.Path(args[args.index('--stats_file') + 1]).write_text('fixture stats\n')
raise SystemExit(status)
PY_VEP
}
"$WORKER_FUNCTION" "$@"
'''


class VEPChromosomeWorkerTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory(prefix="vep_chromosome_worker_",
            dir=os.environ.get("PRIVA_TEST_TMPDIR", os.environ.get("TMPDIR")))
        self.addCleanup(self.temp.cleanup)
        self.folder = Path(self.temp.name).resolve()
        self.cache = self.folder / "cache"
        self.input = self.cache / "clinvar/clinvar.hg38.vep.vcf.gz"
        write_vcf(self.input)
        self.snv = self.cache / "spliceai/snv.vcf.gz"
        self.indel = self.cache / "spliceai/indel.vcf.gz"
        write_vcf(self.snv)
        write_vcf(self.indel)
        self.config = self.folder / "config.yaml"
        self.output = self.folder / "chr1.raw.vcf"
        self.stats = self.folder / "chr1.stats.html"
        self.trace = self.folder / "argv.json"
        self.events = self.folder / "events.txt"
        self.settings = {
            "assembly": "hg38", "vep_installed_version": "116", "threads": 22,
            "cache_root": str(self.cache), "clinvar_vcf_dir": str(self.input.parent),
            "clinvar_vcf": str(self.input), "ref_genome": str(self.cache / "reference.fa"),
            "vep_cache_dir": str(self.cache / "vep"), "vep_plugins_dir": str(self.cache / "plugins"),
            "vep_plugins_cachedir": str(self.cache / "plugin-cache"),
            "loftee_repo": str(self.cache / "loftee"), "loeuf_prescore": "loeuf.tsv.gz",
            "spliceai_snv_prescore": str(self.snv), "spliceai_indel_prescore": str(self.indel),
            "pext_enabled": False,
        }
        for key in ("utr_annotator_file", "alphamissense_prescore", "primateai_prescore",
                    "conservation_file", "splicevault_prescore", "human_ancestor_fasta",
                    "loftee_conservation_file", "gerp_bigwig", "mavedb_file"):
            self.settings[key] = str(self.cache / key)
        self.save_config()

    def save_config(self):
        self.config.write_text(yaml.safe_dump(self.settings, sort_keys=False))

    def run_shell(self, function, arguments, *, oracle=False, env=None):
        environment = dict(os.environ, WORKER_SCRIPT=str(ROOT / "scripts/annotation_vcf.sh"),
            WORKER_ORACLE_PATH=str(ORACLE), WORKER_ORACLE=str(int(oracle)), WORKER_FUNCTION=function,
            WORKER_TRACE=str(self.trace), WORKER_EVENTS=str(self.events), SLURM_CPUS_PER_TASK="24")
        if env:
            environment.update(env)
        return subprocess.run(["bash", "-c", SHELL, "vep-worker-fixture", *map(str, arguments)],
            env=environment, capture_output=True, text=True, timeout=30)

    def worker_args(self, forks=8):
        return ["--config", self.config, "--input-vcf", self.input,
                "--output-vcf", self.output, "--stats-file", self.stats, "--forks", forks]

    def run_worker(self, forks=8, **kwargs):
        return self.run_shell("anno_VEP_cache_worker", self.worker_args(forks), **kwargs)

    def serial(self, *, oracle=False):
        result = self.run_shell("anno_VEP_data", ["--config", self.config,
            "--input_vcf", self.input, "--threads", 22], oracle=oracle)
        self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
        return json.loads(self.trace.read_text()), gzip.open(self.input, "rb").read()

    def reset_worker(self):
        for path in (self.output, self.stats, self.trace, self.events):
            path.unlink(missing_ok=True)

    @staticmethod
    def science_args(argv):
        result = list(argv)
        for flag in ("-i", "-o", "--stats_file", "--fork"):
            result[result.index(flag) + 1] = "<execution>"
        return result

    def test_frozen_oracle_identity(self):
        self.assertEqual(hashlib.sha256(ORACLE.read_bytes()).hexdigest(), ORACLE_SHA)
        function = ORACLE.read_text().split("function _priva_anno_VEP_data() {", 1)[1]
        function = "function _priva_anno_VEP_data() {" + function[:-1]
        self.assertEqual(hashlib.sha256(function.encode()).hexdigest(), FUNCTION_SHA)
        self.assertIn("5be1b9469c06d5ea519c9f66bed2ed286e092c52150c4dc442d979d9769c7f90", ORACLE.read_text())

    def test_serial_22_argv_and_record_output_equal_frozen_v7(self):
        before, old_records = self.serial(oracle=True)
        after, new_records = self.serial()
        self.assertEqual(after, before)
        self.assertEqual(new_records.split(b"#CHROM", 1)[1], old_records.split(b"#CHROM", 1)[1])
        self.assertEqual(after[after.index("--fork") + 1], "22")

    def test_worker_22_control_and_8_candidate_preserve_all_scientific_argv(self):
        original, _ = self.serial(oracle=True)
        for forks in (22, 8):
            with self.subTest(forks=forks):
                self.reset_worker()
                result = self.run_worker(forks)
                self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
                actual = json.loads(self.trace.read_text())
                self.assertEqual(self.science_args(actual), self.science_args(original))
                self.assertEqual(actual[actual.index("--fork") + 1], str(forks))
                self.assertEqual(actual[actual.index("--buffer_size") + 1], "10000")
                self.assertEqual(self.output.read_bytes(), gzip.open(self.input, "rb").read())

    def test_grch37_lof_and_pext_order_equal_frozen_serial(self):
        self.settings["assembly"] = "hg19"
        self.settings["pext_enabled"] = True
        mean = self.folder / "mean.bw"
        mean.write_bytes(b"fixture")
        tissue_root = self.folder / "expression"
        tissues = tissue_root / "tissues"
        tissues.mkdir(parents=True)
        for name in ("Brain-Cortex", "Liver"):
            (tissues / (name + ".bw")).write_bytes(b"fixture")
        self.settings.update(pext_mean_bw=str(mean), pext_tissue_bw_dir=str(tissue_root),
                             pext_tissues="Brain-Cortex,Liver")
        self.save_config()
        original, _ = self.serial(oracle=True)
        self.reset_worker()
        result = self.run_worker()
        self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
        actual = json.loads(self.trace.read_text())
        self.assertEqual(self.science_args(actual), self.science_args(original))
        self.assertFalse(any(item.startswith("MaveDB,") for item in actual))
        lof = next(item for item in actual if item.startswith("LoF,"))
        self.assertIn(",conservation_file:", lof)
        self.assertNotIn(",gerp_bigwig:", lof)
        custom = [actual[i + 1].split(",")[1] for i, value in enumerate(actual[:-1]) if value == "--custom"]
        self.assertEqual(custom, ["PEXT_MEAN", "PEXT_Brain_Cortex", "PEXT_Liver"])

    def test_raw_worker_has_no_skip_sort_index_cleanup_or_global_receipt(self):
        original = self.input.read_bytes()
        config = self.config.read_bytes()
        result = self.run_worker(env={"WORKER_EXISTING_CSQ": "1", "SLURM_CPUS_PER_TASK": "10"})
        self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
        self.assertEqual(self.events.read_text().splitlines(), ["vep"])
        self.assertEqual(self.input.read_bytes(), original)
        self.assertEqual(self.config.read_bytes(), config)
        self.assertTrue(self.stats.is_file())
        for pending in (False, True):
            self.assertFalse(refresh.candidate_receipt_path(self.input, pending=pending).exists())

    def test_failure_status_and_partial_files_retained_without_publication(self):
        original = self.input.read_bytes()
        result = self.run_worker(env={"WORKER_VEP_STATUS": "42"})
        self.assertEqual(result.returncode, 42, result.stdout + result.stderr)
        self.assertEqual(self.events.read_text().splitlines(), ["vep"])
        self.assertEqual(self.output.read_bytes(), gzip.open(self.input, "rb").read()[:100])
        self.assertTrue(self.stats.is_file())
        self.assertEqual(self.input.read_bytes(), original)
        retry = self.run_worker()
        self.assertNotEqual(retry.returncode, 0)
        self.assertEqual(self.events.read_text().splitlines(), ["vep"])

    def test_rejects_missing_duplicate_and_scientific_options_before_vep(self):
        arguments = self.worker_args()
        for invalid in (arguments[:-2], arguments + ["--forks", "9"],
                        arguments + ["--cadd-profile", "unqualified"],
                        arguments + ["--spliceai_snv_prescore", self.snv]):
            with self.subTest(arguments=invalid):
                result = self.run_shell("anno_VEP_cache_worker", invalid)
                self.assertNotEqual(result.returncode, 0)
                self.assertFalse(self.trace.exists())
                self.assertFalse(self.output.exists())

    def test_forks_and_slurm_parent_cpu_budget_validated_before_vep(self):
        for forks, cpus in ((0, "10"), (-1, "10"), ("9.5", "10"), ("x", "10"),
                            ("09", "10"), (9, "10"), (8, "9"), (9, ""), (9, "not-a-number"),
                            ("999999999999999999999999999999", "10")):
            with self.subTest(forks=forks, cpus=cpus):
                result = self.run_worker(forks, env={"SLURM_CPUS_PER_TASK": cpus})
                self.assertNotEqual(result.returncode, 0)
                self.assertFalse(self.trace.exists())
                self.assertFalse(self.output.exists())

    def test_existing_output_or_stats_never_overwritten(self):
        for path in (self.output, self.stats):
            with self.subTest(path=path):
                path.write_bytes(b"retain existing")
                result = self.run_worker()
                self.assertNotEqual(result.returncode, 0)
                self.assertEqual(path.read_bytes(), b"retain existing")
                self.assertFalse(self.trace.exists())
                path.unlink()

    def test_dangling_output_and_parent_alias_rejected(self):
        self.output.symlink_to(self.folder / "missing-target")
        result = self.run_worker()
        self.assertNotEqual(result.returncode, 0)
        self.assertTrue(self.output.is_symlink())
        self.assertFalse(self.trace.exists())
        self.output.unlink()
        alias = self.folder / "alias"
        alias.symlink_to(self.folder, target_is_directory=True)
        self.output = alias / "chr1.raw.vcf"
        result = self.run_worker()
        self.assertNotEqual(result.returncode, 0)
        self.assertFalse(self.trace.exists())

    def test_relative_missing_parent_same_path_and_compressed_output_rejected(self):
        for output in (Path("relative.vcf"), self.folder / "absent/out.vcf",
                       self.stats, self.folder / "wrong.vcf.gz"):
            with self.subTest(output=output):
                self.output = output
                result = self.run_worker()
                self.assertNotEqual(result.returncode, 0)
                self.assertFalse(self.trace.exists())

    def test_worker_explicit_input_output_stats_can_contain_spaces(self):
        source = self.input.parent / "chr1 source.vcf.gz"
        source.write_bytes(self.input.read_bytes())
        self.input = source
        self.output = self.folder / "chr1 output.vcf"
        self.stats = self.folder / "chr1 stats.html"
        result = self.run_worker()
        self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
        self.assertEqual(self.output.read_bytes(), gzip.open(source, "rb").read())
        actual = json.loads(self.trace.read_text())
        self.assertEqual(actual[actual.index("-i") + 1], str(source))
        self.assertEqual(actual[actual.index("-o") + 1], str(self.output))
        self.assertEqual(actual[actual.index("--stats_file") + 1], str(self.stats))

    def test_candidate_phase_owns_only_real_whole_input_receipts(self):
        arguments = ["--config", self.config, "--input-vcf", self.input]
        begin = self.run_shell("anno_VEP_cache_candidate_phase", arguments + ["--phase", "begin"])
        self.assertEqual(begin.returncode, 0, begin.stdout + begin.stderr)
        pending = refresh.candidate_receipt_path(self.input, pending=True)
        self.assertTrue(pending.is_file())
        self.assertFalse(self.trace.exists())
        write_vcf(self.input, annotated=True)
        complete = self.run_shell("anno_VEP_cache_candidate_phase", arguments + ["--phase", "complete"])
        self.assertEqual(complete.returncode, 0, complete.stdout + complete.stderr)
        receipt = json.loads(refresh.candidate_receipt_path(self.input).read_text())
        self.assertFalse(receipt["scientific_qualified"])
        self.assertEqual(receipt["status"], "prepared_annotation_only")
        self.assertFalse(pending.exists())
        self.assertFalse(self.trace.exists())

    def test_phase_rejects_worker_input_missing_begin_and_execution_override(self):
        arguments = ["--config", self.config, "--input-vcf", self.input]
        for suffix in (["--phase", "complete"], ["--phase", "begin", "--forks", "9"],
                       ["--phase", "wrong"], ["--phase", "begin", "--phase", "begin"]):
            with self.subTest(suffix=suffix):
                result = self.run_shell("anno_VEP_cache_candidate_phase", arguments + suffix)
                self.assertNotEqual(result.returncode, 0)
                self.assertFalse(refresh.candidate_receipt_path(self.input).exists())
                self.assertFalse(refresh.candidate_receipt_path(self.input, pending=True).exists())
        shard = self.folder / "chr1.vcf.gz"
        write_vcf(shard)
        result = self.run_shell("anno_VEP_cache_candidate_phase", ["--config", self.config,
            "--input-vcf", shard, "--phase", "begin"])
        self.assertNotEqual(result.returncode, 0)
        self.assertFalse(refresh.candidate_receipt_path(shard, pending=True).exists())

    def test_environment_cannot_turn_private_runtime_into_worker(self):
        # A direct private call retains the original guard even when callers
        # export purpose/output values; this fixture deliberately uses v2.
        write_vcf(self.snv, version2=True)
        write_vcf(self.indel, version2=True)
        result = self.run_shell("_priva_anno_VEP_data", ["--config", self.config,
            "--input_vcf", self.input, "--threads", 9], env={
                "priva_vep_purpose": "cache_worker", "priva_worker_output": str(self.output),
                "priva_worker_stats": str(self.stats)})
        self.assertNotEqual(result.returncode, 0)
        self.assertIn("not qualified", result.stderr)
        self.assertFalse(self.trace.exists())


if __name__ == "__main__":
    unittest.main()
