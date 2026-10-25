"""Real tiny AF split equivalence and the split-to-annotation failure barrier.

python tests/test_gnomad_parallel_splits.py test

Use four allocated CPUs to exercise two unchanged view | sort pipelines. The
retained serial fixture is the scientific oracle; no production cache is read.
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
import tempfile
import unittest


ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "scripts/annotation_vcf.sh"
SERIAL = ROOT / "tests/fixtures/gnomad_serial_e8_0029d817.sh"
SERIAL_SOURCE_SHA256 = "0029d817e6b43792eed944b763094d8ff4acff4bb9872407c904a072dbdc5743"
SERIAL_FUNCTION_SHA256 = "2fd4717defe5ed472ad4d27b698c56c6ab721711bcca4cfaf1f08ac870c959cb"
SPEC = importlib.util.spec_from_file_location(
    "gnomad_split_fixture_inputs", Path(__file__).with_name("test_gnomad_annotation_failures.py"))
FIXTURES = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(FIXTURES)
TOOLS = FIXTURES.TOOLS
COMMAND_HEADER = re.compile(r"^##bcftools_(?:view|sort|annotate|concat)Command=")


SHELL = r'''
source "$1"
[[ -z $TEST_SERIAL_REFERENCE ]] || source "$TEST_SERIAL_REFERENCE"
randomID() { printf 'fixedtag\n'; }
bcftools() {
    if [[ $1 == query && $* == *.chr2.vcf.gz* && $TEST_MODE == split-validation ]]; then
        return 46
    fi
    if [[ $1 == view && $2 == -r ]]; then
        printf 'view-start\t%s\n' "$3" >> "$TEST_EVENTS"
        : > "$TEST_EVENTS.$3.started"
        if [[ $TEST_MODE == success ]]; then
            sleep 0.25
        elif [[ $3 == chr10 && ( $TEST_MODE == split-failure || $TEST_MODE == split-interrupted ) ]]; then
            sleep 0.5
        elif [[ $3 == chr2 && ( $TEST_MODE == split-failure || $TEST_MODE == split-interrupted ) ]]; then
            # Ensure a sibling is genuinely running when this worker fails.
            for ((attempt=0; attempt<200; attempt++)); do
                [[ ! -f $TEST_EVENTS.chr10.started ]] || break
                sleep 0.01
            done
            [[ -f $TEST_EVENTS.chr10.started ]] || return 96
            "$REAL_BCFTOOLS" "$@" | "$REAL_BCFTOOLS" view -i 'ID="chr2_other_allele"' -Ou || return $?
            [[ $TEST_MODE != split-interrupted ]] || kill -TERM "$BASHPID"
            return 42
        fi
        "$REAL_BCFTOOLS" "$@"
        local status=$?
        printf 'view-end\t%s\t%s\n' "$3" "$status" >> "$TEST_EVENTS"
        return "$status"
    fi
    if [[ $1 == annotate ]]; then
        printf 'annotation-start\n' >> "$TEST_EVENTS"
    elif [[ $1 == concat ]]; then
        printf 'concat-start\n' >> "$TEST_EVENTS"
        printf '%s\n' "$@" > "$TEST_CONCAT_ARGS"
    fi
    "$REAL_BCFTOOLS" "$@"
}
parallel() {
    if [[ $* == *anno_agg_gnomAD_split* ]]; then
        printf '%s\n' "$@" > "$TEST_EVENTS.split-command"
        [[ $TEST_MODE != budget-only ]] || return 97
    fi
    "$REAL_PARALLEL" "$@"
}
mv() {
    local source=$1 destination
    [[ $source != -T ]] || source=$2
    for destination in "$@"; do :; done
    if [[ $TEST_MODE == split-publish-index && $source == */.gnomad-split.*/*.chr2.vcf.gz.tbi &&
          $destination == *.chr2.vcf.gz.tbi ]]; then
        return 54
    fi
    command mv "$@"
}
export -f bcftools mv
anno_agg_gnomAD_data "$2" "$TEST_THREADS" hg38 "$3"
status=$?
printf 'parent-exit\t%s\n' "$status" >> "$TEST_EVENTS"
exit "$status"
'''


def physical_pair(path):
    fields = ("st_dev", "st_ino", "st_mode", "st_size", "st_mtime_ns", "st_ctime_ns", "st_nlink")
    return [{"fields": {key: getattr(selected.stat(), key) for key in fields},
             "sha256": hashlib.sha256(selected.read_bytes()).hexdigest()}
            for selected in (path, Path(str(path) + ".tbi"))]


@unittest.skipUnless(all(TOOLS.values()), "real bcftools, htslib and GNU Parallel required")
class ParallelSplitTests(unittest.TestCase):
    info_fields = staticmethod(FIXTURES.AnnotationFailures.info_fields)
    bgzf = FIXTURES.AnnotationFailures.bgzf

    def fixture(self, label, *, stale=False):
        folder, path, gnomad, source, tables = FIXTURES.AnnotationFailures.fixture(self, label, stale=stale)
        # Equal-position alleles, an exact duplicate, an indel and explicit
        # missing values must survive the unchanged default sort/annotation.
        rows = []
        for line in source.splitlines():
            if line.startswith("#CHROM"):
                rows.append('##INFO=<ID=SCORE,Number=1,Type=Float,Description="Missing fixture score">')
            rows.append(line)
            if line.startswith("chr2\t10\t"):
                rows.append(line)
                rows.append("chr2\t10\tchr2_other_allele\tA\tG\t.\tPASS\tORIG=unmatched_allele;SCORE=.")
            elif line.startswith("chr2\t11\t"):
                rows.append("chr2\t11\tchr2_insertion\tA\tAT\t.\tPASS\tORIG=unmatched_insertion")
        source = "\n".join(rows) + "\n"
        path.unlink(); Path(str(path) + ".tbi").unlink()
        self.bgzf(path, source)
        return folder, path, gnomad, source, tables

    def setUp(self):
        temporary = tempfile.TemporaryDirectory(
            prefix="priva_parallel_split_", dir=os.environ.get("PRIVA_TEST_TMPDIR", os.environ.get("TMPDIR")))
        self.addCleanup(temporary.cleanup)
        self.folder = Path(temporary.name)
        self.env = dict(os.environ, LC_ALL="C", SHELL=shutil.which("bash"),
                        REAL_BCFTOOLS=TOOLS["bcftools"], REAL_PARALLEL=TOOLS["parallel"])

    def evidence(self, value):
        parent = os.environ.get("PRIVA_GNOMAD_SPLIT_EVIDENCE_DIR")
        if parent:
            path = Path(parent) / (self._testMethodName + ".json")
            path.parent.mkdir(parents=True, exist_ok=True)
            path.write_text(json.dumps(value, indent=2) + "\n")

    def run_annotation(self, path, gnomad, *, serial=False, mode="success", threads=4, cpus="4"):
        env = dict(self.env, TEST_MODE=mode, TEST_THREADS=str(threads),
                   TEST_SERIAL_REFERENCE=str(SERIAL) if serial else "",
                   TEST_EVENTS=str(path.parent / "events.tsv"),
                   TEST_CONCAT_ARGS=str(path.parent / "concat.args"))
        if cpus is None:
            env.pop("SLURM_CPUS_PER_TASK", None)
        else:
            env["SLURM_CPUS_PER_TASK"] = cpus
        return subprocess.run(["bash", "-c", SHELL, "test", str(SCRIPT), str(path), str(gnomad)],
                              env=env, text=True, capture_output=True, timeout=45)

    def snapshot(self, path):
        header = subprocess.check_output([TOOLS["bcftools"], "view", "--no-version", "-h", str(path)], text=True)
        body = self.body(path)
        contigs = subprocess.check_output([TOOLS["tabix"], "-l", str(path)], text=True).splitlines()
        indexed = {chrom: subprocess.check_output([TOOLS["tabix"], str(path), chrom], text=True)
                   for chrom in contigs}
        return {"scientific_header": [line for line in header.splitlines() if not COMMAND_HEADER.match(line)],
                "command_path_provenance": [line for line in header.splitlines() if COMMAND_HEADER.match(line)],
                "body": body, "index_contigs": contigs, "index_queries": indexed}

    def body(self, path):
        return subprocess.check_output([TOOLS["bcftools"], "view", "--no-version", "-H", str(path)], text=True)

    def test_real_parallel_matches_serial_records_headers_infos_missingness_order_and_index_queries(self):
        self.assertEqual(FIXTURES.sha256(FIXTURES.function_text(SERIAL, "anno_agg_gnomAD_data").encode()),
                         SERIAL_FUNCTION_SHA256)
        outputs, commands, transcripts, raw_splits = {}, {}, {}, {}
        for name in ("serial", "parallel"):
            folder, path, gnomad, _, _ = self.fixture(name, stale=True)
            result = self.run_annotation(path, gnomad, serial=name == "serial")
            self.assertEqual(result.returncode, 0, result.stderr)
            outputs[name] = self.snapshot(path)
            raw_splits[name] = {chrom: self.snapshot(folder / f"input.{chrom}.vcf.gz")
                                for chrom in ("chr2", "chr10", "chr1")}
            commands[name] = (folder / "concat.args").read_text().splitlines()
            transcripts[name] = {"stdout": result.stdout, "stderr": result.stderr}
            if name == "parallel":
                split_command = (folder / "events.tsv.split-command").read_text().splitlines()
                job_log = Path(split_command[split_command.index("--joblog") + 1])
                jobs = [line.split("\t") for line in job_log.read_text().splitlines()[1:]]
                intervals = [(float(row[2]), float(row[2]) + float(row[3])) for row in jobs]
                self.assertEqual(len(jobs), 3)
                self.assertTrue(all(row[6:8] == ["0", "0"] for row in jobs))
                self.assertTrue(any(max(a[0], b[0]) < min(a[1], b[1])
                                    for index, a in enumerate(intervals) for b in intervals[index + 1:]),
                                "The native fixture must execute overlapping split workers")
                active = peak = 0
                for _, change in sorted((stamp, change) for start, end in intervals
                                         for stamp, change in ((start, 1), (end, -1))):
                    active += change
                    peak = max(peak, active)
                self.assertEqual(peak, 2)
                split_measurement = {"native_step": os.environ.get("SLURM_STEP_ID"),
                    "allocated_cpus": os.environ.get("SLURM_CPUS_PER_TASK"),
                    "requested_workers": 2, "maximum_overlapping_workers": peak,
                    "worker_intervals_epoch_seconds": intervals, "raw_joblog_rows": jobs}
                self.assertEqual(split_command[split_command.index("-j") + 1], "2")
                self.assertIn("soon,fail=1", split_command)
        for key in ("scientific_header", "body", "index_contigs", "index_queries"):
            self.assertEqual(outputs["parallel"][key], outputs["serial"][key], key)
            for chrom in ("chr2", "chr10", "chr1"):
                self.assertEqual(raw_splits["parallel"][chrom][key], raw_splits["serial"][chrom][key], (chrom, key))
        normalized_order = {name: [Path(value).name for value in values[2:]] for name, values in commands.items()}
        self.assertEqual(normalized_order["parallel"], normalized_order["serial"])
        self.assertEqual(normalized_order["parallel"], [
            "input.chr1.gnomAD.vcf.gz", "input.chr10.gnomAD.vcf.gz", "input.chr2.gnomAD.vcf.gz"])
        body = outputs["parallel"]["body"].splitlines()
        self.assertEqual(len(body), 9)
        self.assertEqual(sum(line.startswith("chr2\t10\tchr2_10\tA\tC\t") for line in body), 2)
        self.assertIn("chr2\t10\tchr2_other_allele\tA\tG\t.\tPASS\tORIG=unmatched_allele;SCORE=.", body)
        self.assertIn("chr2\t11\tchr2_insertion\tA\tAT\t.\tPASS\tORIG=unmatched_insertion", body)
        self.assertNotIn("unsupported", outputs["parallel"]["body"])
        unmatched = next(line for line in body if line.startswith("chr2\t11\tchr2_11\tA\tC\t"))
        self.assertTrue(unmatched.endswith("\tORIG=keep"), unmatched)
        self.evidence({"outputs": outputs, "raw_splits": raw_splits,
                       "serial_source_sha256": SERIAL_SOURCE_SHA256,
                       "serial_function_sha256": SERIAL_FUNCTION_SHA256,
                       "actual_parallel_measurement": split_measurement,
                       "concat_orders": normalized_order, "transcripts": transcripts,
                       "comparison_exclusions": "Only bcftools view/sort/annotate/concat Command metadata, including command paths/date; all version and scientific header lines retained",
                       "transient_bgzf_and_tbi_byte_identity_claimed": False})

    def test_failed_or_interrupted_split_waits_for_sibling_and_blocks_annotation_and_publication(self):
        results = {}
        for mode, expected in (("split-failure", 42), ("split-interrupted", 143)):
            with self.subTest(mode=mode):
                folder, path, gnomad, _, _ = self.fixture(mode)
                before = physical_pair(path)
                old = folder / ".gnomad-split.retained-old"
                old.mkdir(); sentinel = old / "input.chr2.vcf.gz"; sentinel.write_bytes(b"retain old hidden stage")
                result = self.run_annotation(path, gnomad, mode=mode)
                events = (folder / "events.tsv").read_text().splitlines()
                self.assertEqual(result.returncode, expected, result.stderr)
                self.assertEqual(physical_pair(path), before)
                self.assertNotIn("annotation-start", events)
                self.assertNotIn("concat-start", events)
                self.assertFalse((folder / "input.chr2.vcf.gz").exists())
                self.assertFalse(Path(str(path) + ".gnomAD.publish.pending").exists())
                self.assertLess(events.index("view-end\tchr10\t0"), events.index(f"parent-exit\t{expected}"))
                self.assertEqual(sentinel.read_bytes(), b"retain old hidden stage")
                partials = [p for p in folder.glob(".gnomad-split.*/input.chr2.vcf.gz") if p != sentinel]
                self.assertEqual(len(partials), 1)
                self.assertEqual(len(self.body(partials[0]).splitlines()), 1)
                retried = self.run_annotation(path, gnomad)
                self.assertEqual(retried.returncode, 0, retried.stderr)
                self.assertEqual(len(self.snapshot(path)["body"].splitlines()), 9)
                self.assertTrue(partials[0].is_file())
                self.assertEqual(sentinel.read_bytes(), b"retain old hidden stage")
                results[mode] = {"failed_status": result.returncode, "retry_status": retried.returncode,
                                 "events": events, "stderr": result.stderr, "original_pair_unchanged": True,
                                 "partial_stage_retained": True, "annotation_and_concat_not_started": True}
        self.evidence(results)

    def test_split_validation_or_index_publication_failure_blocks_annotation(self):
        evidence = {}
        for mode, expected in (("split-validation", 10), ("split-publish-index", 54)):
            with self.subTest(mode=mode):
                folder, path, gnomad, _, _ = self.fixture(mode)
                before = physical_pair(path)
                result = self.run_annotation(path, gnomad, mode=mode)
                self.assertEqual(result.returncode, expected, result.stderr)
                self.assertEqual(physical_pair(path), before)
                events = (folder / "events.tsv").read_text().splitlines()
                self.assertNotIn("annotation-start", events)
                self.assertNotIn("concat-start", events)
                split = folder / "input.chr2.vcf.gz"
                self.assertEqual(split.exists(), mode == "split-publish-index")
                self.assertFalse(Path(str(split) + ".tbi").exists())
                self.assertTrue(list(folder.glob(".gnomad-split.*/*chr2.vcf.gz*")))
                retried = self.run_annotation(path, gnomad)
                self.assertEqual(retried.returncode, 0, retried.stderr)
                self.assertEqual(len(self.snapshot(path)["body"].splitlines()), 9)
                evidence[mode] = {"failed_status": result.returncode, "retry_status": retried.returncode,
                                  "original_pair_unchanged_after_failure": True,
                                  "annotation_and_concat_not_started": True, "stderr": result.stderr}
        self.evidence(evidence)

    def test_existing_complete_split_and_index_are_reused_unchanged(self):
        folder, path, gnomad, _, _ = self.fixture("reuse")
        prepared = subprocess.run(["bash", str(SCRIPT), "anno_agg_gnomAD_split",
            "--input-vcf", str(path), "--chromosome", "chr2"], env=self.env,
            text=True, capture_output=True, timeout=20)
        self.assertEqual(prepared.returncode, 0, prepared.stderr)
        split = folder / "input.chr2.vcf.gz"; before = physical_pair(split)
        result = self.run_annotation(path, gnomad)
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertEqual(physical_pair(split), before)
        self.assertNotIn("view-start\tchr2", (folder / "events.tsv").read_text().splitlines())
        self.evidence({"returncode": result.returncode, "split_pair_unchanged": True,
                       "no_repeated_view_sort": True, "stderr": result.stderr})

    def test_config_and_slurm_budget_bound_workers_with_serial_minimum(self):
        cases = ((1, None, 1), (3, None, 1), (4, None, 2), (22, None, 11),
                 (22, "4", 2), (4, "22", 2), (22, "3", 1), (22, "1", 1))
        evidence = []
        for index, (threads, cpus, expected) in enumerate(cases):
            with self.subTest(threads=threads, cpus=cpus):
                folder, path, gnomad, _, _ = self.fixture("budget-" + str(index))
                result = self.run_annotation(path, gnomad, mode="budget-only", threads=threads, cpus=cpus)
                self.assertEqual(result.returncode, 97, result.stderr)
                command = (folder / "events.tsv.split-command").read_text().splitlines()
                self.assertEqual(command[command.index("-j") + 1], str(expected))
                self.assertNotIn("annotation-start", (folder / "events.tsv").read_text().splitlines())
                evidence.append({"config_threads": threads, "slurm_cpus": cpus, "workers": expected,
                                 "mock_dispatch_no_workers_executed": True})
        self.evidence(evidence)


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    commands = parser.add_subparsers(dest="procedure", required=True)
    run = commands.add_parser("test")
    run.add_argument("--verbosity", type=int, default=2)
    args = parser.parse_args(argv)
    suite = unittest.defaultTestLoader.loadTestsFromTestCase(ParallelSplitTests)
    result = unittest.TextTestRunner(verbosity=args.verbosity).run(suite)
    return 0 if result.wasSuccessful() else 1


if __name__ == "__main__":
    raise SystemExit(main())
