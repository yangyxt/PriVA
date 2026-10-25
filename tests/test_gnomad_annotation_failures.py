"""Fail closed on gnomAD annotation errors; retain the original successful result."""
import hashlib
import json
import os
from pathlib import Path
import shutil
import subprocess
import tempfile
import unittest

ROOT = Path(__file__).resolve().parents[1]
ANNOTATION = ROOT / "scripts/annotation_vcf.sh"
COMMON = ROOT / "scripts/common_bash_utils.sh"
REFERENCE = Path(os.environ.get(
    "PRIVA_GNOMAD_ANNOTATION_REFERENCE_SCRIPTS",
    ROOT / "tests/fixtures/installer-054b494/scripts"))
TOOLS = {name: shutil.which(name) for name in ("bcftools", "bgzip", "tabix", "parallel")}
HEADER = "Seq\tHost\tStarttime\tJobRuntime\tSend\tReceive\tExitval\tSignal\tCommand\n"
ORIGINAL_FUNCTION_HASHES = {
    "annotation_vcf.sh": ("anno_agg_gnomAD_data", "3393babd3734c22f65b095004fa853f3735e24ac68b3c300d92c35fb668ede65"),
    "common_bash_utils.sh": ("check_parallel_joblog", "66168c157113163b360873ef4974969d5ca956e7515c8516d24a11fdcaed8d12"),
}


def sha256(data):
    return hashlib.sha256(data).hexdigest()


def function_text(path, name):
    text = path.read_text()
    start = text.index("function " + name)
    return text[start:text.index("\n}\n", start) + 3]


def record_evidence(label, value):
    target = os.environ.get("PRIVA_GNOMAD_EVIDENCE_DIR")
    if target:
        path = Path(target) / (label + ".json")
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(json.dumps(value, indent=2) + "\n")


class JoblogChecks(unittest.TestCase):
    def setUp(self):
        temporary = tempfile.TemporaryDirectory(
            prefix="priva_joblog_", dir=os.environ.get("PRIVA_TEST_TMPDIR", os.environ.get("TMPDIR")))
        self.addCleanup(temporary.cleanup)
        self.folder = Path(temporary.name)

    def check(self, content, expected=1, status=0, legacy=False):
        path = self.folder / "jobs.log"
        if content is not None:
            path.write_text(content)
        call = ('check_parallel_joblog "$2"' if legacy else
                'check_parallel_joblog --job-log "$2" --expected-jobs "$3" --exit-status "$4"')
        setup = 'source "$1"\n'
        if legacy:
            setup += '(exit "$4")\n'
        result = subprocess.run(
            ["bash", "-c", setup + call, "test", str(COMMON), str(path), str(expected), str(status)],
            text=True, capture_output=True, timeout=15)
        record_evidence(self._testMethodName, {
            "returncode": result.returncode, "stdout": result.stdout, "stderr": result.stderr,
            "joblog": content, "expected_jobs": expected, "incoming_status": status,
            "legacy_interface": legacy,
        })
        return result

    def test_valid_joblog_and_legacy_single_path_interface(self):
        content = HEADER + "1\t:\t1.1\t0.1\t0\t0\t0\t0\techo ready\n"
        self.assertEqual(self.check(content).returncode, 0)
        self.assertEqual(self.check(content, legacy=True).returncode, 0)

    def test_prior_parallel_status_is_retained_even_with_missing_log(self):
        for status in (23, 137):
            with self.subTest(status=status):
                self.assertEqual(self.check(None, status=status).returncode, status)
                self.assertEqual(self.check(None, status=status, legacy=True).returncode, status)

    def test_256_failed_jobs_cannot_wrap_to_success(self):
        rows = "".join(f"{n}\t:\t1\t0\t0\t0\t1\t0\tfalse\n" for n in range(1, 257))
        self.assertEqual(self.check(HEADER + rows, expected=256).returncode, 1)

    def test_nonzero_exit_and_signal_are_retained(self):
        self.assertEqual(self.check(HEADER + "1\t:\t1\t0\t0\t0\t17\t0\tfalse\n").returncode, 17)
        self.assertEqual(self.check(HEADER + "1\t:\t1\t0\t0\t0\t0\t9\tsleep 1\n").returncode, 137)

    def test_missing_empty_malformed_and_incomplete_logs_fail(self):
        malformed = [
            None, "", HEADER, "not a GNU Parallel joblog\n",
            HEADER + "1\t:\t1\t0\t0\t0\t0\n",
            HEADER + "1\t:\t1\t0\t0\t0\twrong\t0\techo ready\n",
            HEADER + "1\t:\t1\t0\t0\t0\t0\twrong\techo ready\n",
            HEADER + "1\t:\t1\t0\t0\t0\t256\t0\techo ready\n",
            HEADER + "1\t:\t1\t0\t0\t0\t0\t128\techo ready\n",
            HEADER + "1\t:\t1\t0\t0\t0\t0\t0\techo ready\n" * 2,
            HEADER + "2\t:\t1\t0\t0\t0\t0\t0\techo ready\n",
        ]
        for content in malformed:
            with self.subTest(content=content):
                path = self.folder / "jobs.log"
                if path.exists():
                    path.unlink()
                self.assertNotEqual(self.check(content).returncode, 0)
        content = HEADER + "1\t:\t1\t0\t0\t0\t0\t0\techo ready\n"
        self.assertNotEqual(self.check(content, expected=2).returncode, 0)


SHELL = r'''
source "$1"
randomID() { printf 'fixedtag\n'; }
eval "$(declare -f check_vcf_validity | sed '1s/check_vcf_validity/original_check_vcf_validity/')"
check_vcf_validity() {
    case "$TEST_FAILURE:$1" in
        split-validation:*.chr2.vcf.gz) return 46 ;;
        annotated-validation:*.chr2.gnomAD.vcf.gz) return 47 ;;
        final-validation:*.fixedtag.vcf.gz) return 48 ;;
    esac
    original_check_vcf_validity "$@"
}
bcftools() {
    # Split validation now runs in its own modular worker, which sources the
    # real common validator. Inject its tool failure across that process edge.
    if [[ $1 == query && $* == *.chr2.vcf.gz* && $TEST_FAILURE == split-validation ]]; then
        return 46
    fi
    if [[ $1 == query && $2 == -f && $TEST_FAILURE == query ]]; then
        printf 'chr1\n'
        return 41
    fi
    if [[ $1 == view && $2 == -r && $3 == chr2 && $TEST_FAILURE == split-view ]]; then
        "$REAL_BCFTOOLS" "$@" | "$REAL_BCFTOOLS" view -i "POS=10" -Ou
        return 42
    fi
    if [[ $1 == annotate && $* == *.chr2.gnomAD.vcf.gz* &&
          ( $TEST_FAILURE == annotation-producer || $TEST_FAILURE == annotation-interrupted ) ]]; then
        "$REAL_BCFTOOLS" "$@" || return $?
        local output="" previous="" argument
        for argument in "$@"; do
            [[ $previous != -o ]] || output=$argument
            previous=$argument
        done
        "$REAL_BCFTOOLS" view -i "POS=10" -Oz -o "$output.one" "$output" &&
        mv "$output.one" "$output" || return $?
        [[ $TEST_FAILURE != annotation-interrupted ]] || kill -TERM "$BASHPID"
        return 53
    fi
    if [[ $1 == sort ]]; then
        case "$TEST_FAILURE:$*" in
            split-sort:*.chr2.vcf.gz*) return 43 ;;
            merge-sort:*.fixedtag.vcf.gz*) return 50 ;;
        esac
    fi
    if [[ $1 == concat ]]; then
        printf '%s\n' "$@" > "$TEST_CONCAT_ARGS"
        "$REAL_BCFTOOLS" "$@" || return $?
        [[ $TEST_FAILURE != merge-concat ]] || return 49
        return 0
    fi
    "$REAL_BCFTOOLS" "$@"
}
tabix() {
    case "$TEST_FAILURE:$*" in
        split-tabix:*.chr2.vcf.gz) return 44 ;;
        final-tabix:*.fixedtag.vcf.gz) return 45 ;;
    esac
    case "$TEST_FAILURE:$*" in
        repair-index:*/.gnomad-index.*/*) return 55 ;;
    esac
    "$REAL_TABIX" "$@"
}
parallel() {
    # These existing parallel/joblog faults exercise the annotation barrier;
    # split-worker faults have their own real-pipeline checks below.
    if [[ $* == *anno_agg_gnomAD_split* ]]; then
        "$REAL_PARALLEL" "$@"
        return $?
    fi
    local job_log="" previous="" argument dry=0
    for argument in "$@"; do
        [[ $previous != --joblog ]] || job_log=$argument
        [[ $argument != --dry-run ]] || dry=1
        previous=$argument
    done
    if [[ $dry == 1 ]]; then
        [[ $TEST_FAILURE != parallel-dry ]] || return 51
        "$REAL_PARALLEL" "$@"
        return $?
    fi
    "$REAL_PARALLEL" "$@" || return $?
    case "$TEST_FAILURE" in
        parallel-run) return 52 ;;
        joblog-missing) rm -f "$job_log" ;;
        joblog-empty) : > "$job_log" ;;
        joblog-malformed) printf 'broken\n' > "$job_log" ;;
        joblog-incomplete) head -n 2 "$job_log" > "${job_log}.one"; mv "${job_log}.one" "$job_log" ;;
        joblog-nonzero) awk -F '\t' 'BEGIN{OFS=FS} NR==2{$7=17}{print}' "$job_log" > "${job_log}.changed"; mv "${job_log}.changed" "$job_log" ;;
        joblog-signal) awk -F '\t' 'BEGIN{OFS=FS} NR==2{$8=9}{print}' "$job_log" > "${job_log}.changed"; mv "${job_log}.changed" "$job_log" ;;
    esac
}
mv() {
    local source=$1 destination
    [[ $source != -T ]] || source=$2
    for destination in "$@"; do :; done
    if [[ $TEST_FAILURE == publish-index && $source == *.fixedtag.vcf.gz.tbi &&
          $destination == "$TEST_INPUT.tbi" ]]; then
        return 54
    fi
    if [[ $TEST_FAILURE == chromosome-publish-index &&
          $source == */.gnomad-annotation.*/* && $destination == *.chr2.gnomAD.vcf.gz.tbi ]]; then
        return 54
    fi
    command mv "$@"
}
export -f bcftools tabix mv
anno_agg_gnomAD_data "$2" 1 hg38 "$3"
'''


@unittest.skipUnless(all(TOOLS.values()), "real bcftools, htslib and GNU Parallel required")
class AnnotationFailures(unittest.TestCase):
    def setUp(self):
        temporary = tempfile.TemporaryDirectory(
            prefix="priva_gnomad_annotation_", dir=os.environ.get("PRIVA_TEST_TMPDIR", os.environ.get("TMPDIR")))
        self.addCleanup(temporary.cleanup)
        self.folder = Path(temporary.name)
        self.env = dict(os.environ, LC_ALL="C", SHELL=shutil.which("bash"), TEST_FAILURE="", REAL_BCFTOOLS=TOOLS["bcftools"],
                        REAL_TABIX=TOOLS["tabix"], REAL_PARALLEL=TOOLS["parallel"])

    @staticmethod
    def info_fields():
        fields = ["AC_joint", "AN_joint", "AF_joint", "nhomalt_joint",
                  "AC_grpmax_joint", "AF_grpmax_joint", "AN_grpmax_joint", "nhomalt_grpmax_joint",
                  "AC_joint_XX", "AN_joint_XX", "AC_joint_XY", "AN_joint_XY",
                  "AF_joint_XX", "AF_joint_XY", "nhomalt_joint_XX", "nhomalt_joint_XY"]
        for population in ("nfe", "eas", "afr", "amr", "asj", "fin", "sas", "mid", "remaining"):
            fields += [f"AF_joint_{population}_XX", f"AF_joint_{population}_XY"]
        fields += [f"nhomalt_joint_{p}_XX" for p in
                   ("nfe", "eas", "afr", "amr", "asj", "fin", "sas", "mid", "remaining")]
        return fields

    def bgzf(self, path, text):
        with path.open("wb") as output:
            subprocess.run([TOOLS["bgzip"], "-c"], input=text.encode(), stdout=output, check=True)
        subprocess.run([TOOLS["tabix"], "-p", "vcf", str(path)], check=True)
        return path

    def fixture(self, label, *, stale=False, unsupported_only=False):
        folder = self.folder / label
        folder.mkdir()
        header = "##fileformat=VCFv4.2\n"
        for chrom in ("chr2", "chr10", "chr1", "chr3", "GL000191.1"):
            header += f"##contig=<ID={chrom},length=1000>\n"
        header += '##INFO=<ID=ORIG,Number=1,Type=String,Description="Original fixture marker">\n'
        columns = "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
        rows = ""
        if not unsupported_only:
            for chrom in ("chr2", "chr10", "chr1"):
                for position in (10, 11):
                    rows += f"{chrom}\t{position}\t{chrom}_{position}\tA\tC\t50\tPASS\tORIG=keep\n"
        rows += "GL000191.1\t10\tunsupported\tA\tC\t50\tPASS\tORIG=skip_original_policy\n"
        source = header + columns + rows
        path = self.bgzf(folder / "input.vcf.gz", source)
        fields = self.info_fields()
        annotated_header = header
        for field in fields:
            kind = "Float" if field.startswith("AF_") else "Integer"
            annotated_header += f'##INFO=<ID={field},Number=1,Type={kind},Description="Synthetic {field}">\n'
        tables = {}
        for chrom in ("chr1", "chr2", "chr3", "chr10", "chrX"):
            records = ""
            for position in ((10,) if chrom == "chr2" else (10, 11)):
                info = ";".join(f"{field}={'0.125' if field.startswith('AF_') else '2'}" for field in fields)
                records += f"{chrom}\t{position}\t.\tA\tC\t.\tPASS\t{info}\n"
            text = annotated_header
            if chrom == "chrX":
                text += "##contig=<ID=chrX,length=1000>\n"
            text += columns + records
            tables[chrom] = text
            self.bgzf(folder / f"gnomad.chr{chrom[3:]}.vcf.gz", text)
        if stale:
            self.bgzf(folder / "input.chr3.gnomAD.vcf.gz", tables["chr3"])
        return folder, path, folder / "gnomad.chrX.vcf.gz", source, tables

    def run_annotation(self, path, gnomad, *, failure="", source=ANNOTATION):
        env = dict(self.env, TEST_FAILURE=failure, TEST_INPUT=str(path),
                   TEST_CONCAT_ARGS=str(path.parent / "concat.args"))
        return subprocess.run(["bash", "-c", SHELL, "test", str(source), str(path), str(gnomad)],
                              env=env, text=True, capture_output=True, timeout=75)

    def assert_failure(self, failure, expected):
        folder, path, gnomad, _, _ = self.fixture(failure, stale=True)
        before = path.read_bytes(), Path(str(path) + ".tbi").read_bytes()
        stale = folder / "input.chr3.gnomAD.vcf.gz"
        stale_before = stale.read_bytes()
        result = self.run_annotation(path, gnomad, failure=failure)
        unchanged = before == (path.read_bytes(), Path(str(path) + ".tbi").read_bytes())
        record_evidence(self._testMethodName, {
            "failure": failure, "expected_returncode": expected, "returncode": result.returncode,
            "selected_vcf_and_index_unchanged": unchanged,
            "input_sha256": sha256(before[0]), "after_sha256": sha256(path.read_bytes()),
            "stdout": result.stdout, "stderr": result.stderr,
        })
        self.assertEqual(result.returncode, expected, result.stderr)
        self.assertTrue(unchanged, "selected VCF/index were published despite failure")
        self.assertEqual(stale.read_bytes(), stale_before)

    def test_complete_real_header_reuses_without_entering_annotation(self):
        folder, path, gnomad, _, _ = self.fixture("warm")
        cold = self.run_annotation(path, gnomad)
        self.assertEqual(cold.returncode, 0, cold.stderr)
        before = path.read_bytes(), Path(str(path) + ".tbi").read_bytes()
        # A query fault can only fire after the header readiness branch.
        warm = self.run_annotation(path, gnomad, failure="query")
        self.assertEqual(warm.returncode, 0, warm.stderr)
        self.assertIn("We do not need to add them again", warm.stderr + warm.stdout)
        self.assertEqual(before, (path.read_bytes(), Path(str(path) + ".tbi").read_bytes()))
        old = self.run_annotation(path, gnomad, failure="query", source=REFERENCE / "annotation_vcf.sh")
        self.assertNotIn("We do not need to add them again", old.stderr + old.stdout)

    def test_missing_required_header_field_enters_cold_path_and_propagates_failure(self):
        for field in ("AF_grpmax_joint", "AN_joint_XY", "nhomalt_joint_remaining_XX"):
            with self.subTest(field=field):
                folder, path, gnomad, _, _ = self.fixture(field)
                self.assertEqual(self.run_annotation(path, gnomad).returncode, 0)
                stripped = folder / "stripped.vcf.gz"
                subprocess.run([TOOLS["bcftools"], "annotate", "-x", "INFO/" + field,
                                "-Oz", "-o", str(stripped), str(path)], check=True)
                stripped.replace(path)
                subprocess.run([TOOLS["tabix"], "-f", "-p", "vcf", str(path)], check=True)
                before = path.read_bytes(), Path(str(path) + ".tbi").read_bytes()
                result = self.run_annotation(path, gnomad, failure="query")
                self.assertEqual(result.returncode, 41, result.stderr)
                self.assertEqual(before, (path.read_bytes(), Path(str(path) + ".tbi").read_bytes()))

    def test_real_success_matches_pinned_original_rows_infos_and_order(self):
        for name, (function, checksum) in ORIGINAL_FUNCTION_HASHES.items():
            self.assertEqual(sha256(function_text(REFERENCE / name, function).encode()), checksum)
        outputs = {}
        transcripts = {}
        source_text = tables = None
        for label, script in (("original", REFERENCE / "annotation_vcf.sh"), ("fixed", ANNOTATION)):
            folder, path, gnomad, source_text, tables = self.fixture(label)
            result = self.run_annotation(path, gnomad, source=script)
            self.assertEqual(result.returncode, 0, result.stderr)
            body = subprocess.check_output(
                [TOOLS["bcftools"], "query", "-f", "%CHROM\t%POS\t%ID\t%REF\t%ALT\t%QUAL\t%FILTER\t%INFO\n", str(path)],
                text=True)
            header = subprocess.check_output([TOOLS["bcftools"], "view", "--no-version", "-h", str(path)], text=True)
            header = "\n".join(line for line in header.splitlines() if not line.startswith("##bcftools"))
            concat = [Path(x).name for x in (folder / "concat.args").read_text().splitlines()[2:]]
            outputs[label] = {"body": body, "header": header, "concat_order": concat,
                              "output_sha256": sha256(path.read_bytes())}
            transcripts[label] = {"stdout": result.stdout, "stderr": result.stderr}
        record_evidence(self._testMethodName, {
            "reference_scripts": str(REFERENCE), "fixed_scripts": str(ANNOTATION.parent),
            "reference_function_hashes": ORIGINAL_FUNCTION_HASHES,
            "reference_source_hashes": {n: sha256((REFERENCE / n).read_bytes()) for n in ORIGINAL_FUNCTION_HASHES},
            "fixed_source_hashes": {n: sha256((ANNOTATION.parent / n).read_bytes()) for n in ORIGINAL_FUNCTION_HASHES},
            "input_vcf_text": source_text, "gnomad_vcf_texts": tables,
            "outputs": outputs, "transcripts": transcripts,
            "same_body": outputs["original"]["body"] == outputs["fixed"]["body"],
            "same_header": outputs["original"]["header"] == outputs["fixed"]["header"],
            "same_concat_order": outputs["original"]["concat_order"] == outputs["fixed"]["concat_order"],
        })
        self.assertEqual(outputs["fixed"]["body"], outputs["original"]["body"])
        self.assertEqual(outputs["fixed"]["header"], outputs["original"]["header"])
        self.assertEqual(outputs["fixed"]["concat_order"], outputs["original"]["concat_order"])
        self.assertEqual(outputs["fixed"]["concat_order"], [
            "input.chr1.gnomAD.vcf.gz", "input.chr10.gnomAD.vcf.gz", "input.chr2.gnomAD.vcf.gz"])
        self.assertEqual(len(outputs["fixed"]["body"].splitlines()), 6)
        self.assertNotIn("unsupported", outputs["fixed"]["body"])
        unmatched = next(line for line in outputs["fixed"]["body"].splitlines() if line.startswith("chr2\t11\t"))
        self.assertTrue(unmatched.endswith("\tORIG=keep"))

    def test_failed_and_interrupted_partial_shards_are_not_reused_on_retry(self):
        for failure in ("split-view", "annotation-producer", "annotation-interrupted"):
            with self.subTest(failure=failure):
                folder, path, gnomad, _, _ = self.fixture("retry-" + failure)
                before = path.read_bytes(), Path(str(path) + ".tbi").read_bytes()
                failed = self.run_annotation(path, gnomad, failure=failure)
                self.assertNotEqual(failed.returncode, 0, failed.stderr)
                self.assertEqual(before, (path.read_bytes(), Path(str(path) + ".tbi").read_bytes()))
                name = "input.chr2.vcf.gz" if failure == "split-view" else "input.chr2.gnomAD.vcf.gz"
                self.assertFalse((folder / name).exists(), "failed producer left a reusable canonical shard")
                stages = list(folder.glob(".gnomad-*/*" + name))
                self.assertTrue(stages, failed.stderr)
                staged_rows = subprocess.check_output(
                    [TOOLS["bcftools"], "query", "-f", "%CHROM\t%POS\n", str(stages[0])], text=True)
                self.assertEqual(len(staged_rows.splitlines()), 1, "fixture must leave a plausible partial shard")
                retried = self.run_annotation(path, gnomad)
                self.assertEqual(retried.returncode, 0, retried.stderr)
                body = subprocess.check_output(
                    [TOOLS["bcftools"], "query", "-f", "%CHROM\t%POS\n", str(path)], text=True)
                self.assertEqual(len(body.splitlines()), 6)
                self.assertEqual(sum(line.startswith("chr2\t") for line in body.splitlines()), 2)
                record_evidence(self._testMethodName + "-" + failure, {
                    "failure_returncode": failed.returncode, "retry_returncode": retried.returncode,
                    "failed_stdout": failed.stdout, "failed_stderr": failed.stderr,
                    "retry_stdout": retried.stdout, "retry_stderr": retried.stderr,
                    "partial_staged_rows": staged_rows, "successful_retry_rows": body,
                    "canonical_failed_shard_absent_before_retry": True,
                    "selected_input_unchanged_after_failure": True,
                })

    def test_failed_final_index_rename_keeps_marker_until_staged_repair_succeeds(self):
        folder, path, gnomad, _, _ = self.fixture("publish-index")
        index = Path(str(path) + ".tbi")
        pending = Path(str(path) + ".gnomAD.publish.pending")
        original_index = index.read_bytes()
        failed = self.run_annotation(path, gnomad, failure="publish-index")
        self.assertEqual(failed.returncode, 54, failed.stderr)
        self.assertTrue(pending.is_file())
        self.assertEqual(index.read_bytes(), original_index)
        completed_data = path.read_bytes()
        expected_rows = subprocess.check_output(
            [TOOLS["bcftools"], "query", "-f", "%CHROM\t%POS\t%INFO\n", str(path)], text=True)
        self.assertEqual(len(expected_rows.splitlines()), 6)
        repair_failed = self.run_annotation(path, gnomad, failure="repair-index")
        self.assertEqual(repair_failed.returncode, 55, repair_failed.stderr)
        self.assertEqual(index.read_bytes(), original_index, "failed repair overwrote the previous index")
        self.assertEqual(path.read_bytes(), completed_data)
        self.assertTrue(pending.exists(), "failed repair removed the pending publication marker")
        retried = self.run_annotation(path, gnomad)
        self.assertEqual(retried.returncode, 0, retried.stderr)
        self.assertFalse(pending.exists())
        indexed = subprocess.check_output(
            [TOOLS["bcftools"], "view", "--no-version", "-H", "-r", "chr2:10-11", str(path)], text=True)
        self.assertEqual(len(indexed.splitlines()), 2)
        actual_rows = subprocess.check_output(
            [TOOLS["bcftools"], "query", "-f", "%CHROM\t%POS\t%INFO\n", str(path)], text=True)
        self.assertEqual(actual_rows, expected_rows)
        self.assertNotEqual(index.read_bytes(), original_index)
        record_evidence(self._testMethodName, {
            "failed_publication_returncode": failed.returncode,
            "failed_repair_returncode": repair_failed.returncode, "retry_returncode": retried.returncode,
            "old_index_preserved_through_failed_repair": True,
            "completed_data_preserved_through_failed_repair": True,
            "marker_retained_through_failed_repair": True, "marker_cleared_after_success": True,
            "expected_rows": expected_rows, "actual_rows": actual_rows, "indexed_chr2_rows": indexed,
            "publication_stderr": failed.stderr, "repair_stderr": repair_failed.stderr,
            "retry_stderr": retried.stderr,
        })

    def test_failed_chromosome_index_rename_is_repaired_before_shard_reuse(self):
        folder, path, gnomad, _, _ = self.fixture("chromosome-index")
        before = path.read_bytes(), Path(str(path) + ".tbi").read_bytes()
        shard = folder / "input.chr2.gnomAD.vcf.gz"
        pending = Path(str(shard) + ".gnomAD.publish.pending")
        failed = self.run_annotation(path, gnomad, failure="chromosome-publish-index")
        self.assertNotEqual(failed.returncode, 0, failed.stderr)
        self.assertEqual(before, (path.read_bytes(), Path(str(path) + ".tbi").read_bytes()))
        self.assertTrue(shard.is_file())
        self.assertFalse(Path(str(shard) + ".tbi").exists())
        self.assertTrue(pending.is_file())
        retried = self.run_annotation(path, gnomad)
        self.assertEqual(retried.returncode, 0, retried.stderr)
        self.assertFalse(pending.exists())
        indexed = subprocess.check_output(
            [TOOLS["bcftools"], "view", "--no-version", "-H", "-r", "chr2:10-11", str(path)], text=True)
        self.assertEqual(len(indexed.splitlines()), 2)
        record_evidence(self._testMethodName, {
            "failed_publication_returncode": failed.returncode, "retry_returncode": retried.returncode,
            "selected_input_unchanged_after_shard_failure": True,
            "pending_shard_index_repaired_before_reuse": True, "indexed_chr2_rows": indexed,
            "publication_stderr": failed.stderr, "retry_stderr": retried.stderr,
        })

    def test_linked_or_directory_publication_markers_are_rejected(self):
        for kind in ("symlink", "directory"):
            for selected in ("input.vcf.gz", "input.chr2.gnomAD.vcf.gz"):
                with self.subTest(kind=kind, selected=selected):
                    folder, path, gnomad, _, _ = self.fixture(kind + selected)
                    target = folder / "untouched-marker-target"
                    target.write_bytes(b"retain original marker target")
                    marker = folder / (selected + ".gnomAD.publish.pending")
                    if kind == "symlink":
                        marker.symlink_to(target)
                    else:
                        marker.mkdir()
                    before = path.read_bytes(), Path(str(path) + ".tbi").read_bytes()
                    result = self.run_annotation(path, gnomad)
                    self.assertNotEqual(result.returncode, 0)
                    self.assertEqual(before, (path.read_bytes(), Path(str(path) + ".tbi").read_bytes()))
                    self.assertEqual(target.read_bytes(), b"retain original marker target")

    def test_directory_index_destination_is_not_reused_or_treated_as_publication(self):
        for selected in ("input.vcf.gz", "input.chr2.vcf.gz", "input.chr2.gnomAD.vcf.gz"):
            with self.subTest(selected=selected):
                folder, path, gnomad, _, _ = self.fixture("directory-" + selected)
                destination = folder / (selected + ".tbi")
                if destination.exists():
                    destination.unlink()
                destination.mkdir()
                sentinel = destination / "retain.txt"
                sentinel.write_bytes(b"retain directory contents")
                before = path.read_bytes()
                result = self.run_annotation(path, gnomad)
                self.assertNotEqual(result.returncode, 0)
                self.assertEqual(path.read_bytes(), before)
                self.assertEqual(list(destination.iterdir()), [sentinel])
                self.assertEqual(sentinel.read_bytes(), b"retain directory contents")
                record_evidence(self._testMethodName + "-" + selected, {
                    "index_destination": selected + ".tbi", "returncode": result.returncode,
                    "selected_data_unchanged": True, "directory_contents_unchanged": True,
                    "stdout": result.stdout, "stderr": result.stderr,
                })

    def test_valid_symlink_to_file_index_remains_usable_without_replacement(self):
        folder, path, gnomad, _, _ = self.fixture("index-symlink")
        index = Path(str(path) + ".tbi")
        original = folder / "original.tbi"
        index.replace(original)
        index.symlink_to(original)
        before = path.read_bytes(), original.read_bytes()
        result = subprocess.run(
            ["bash", str(ANNOTATION), "anno_agg_gnomAD_index", "--input-vcf", str(path)],
            env=self.env, text=True, capture_output=True, timeout=20)
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertTrue(index.is_symlink())
        self.assertEqual(before, (path.read_bytes(), original.read_bytes()))
        record_evidence(self._testMethodName, {
            "returncode": result.returncode, "valid_index_link_preserved": True,
            "selected_data_and_original_index_unchanged": True,
            "stdout": result.stdout, "stderr": result.stderr,
        })

    def test_unrelated_stale_shard_is_neither_published_nor_deleted(self):
        folder, path, gnomad, _, _ = self.fixture("stale", stale=True)
        stale = folder / "input.chr3.gnomAD.vcf.gz"
        original = stale.read_bytes(), Path(str(stale) + ".tbi").read_bytes()
        result = self.run_annotation(path, gnomad)
        self.assertEqual(result.returncode, 0, result.stderr)
        body = subprocess.check_output([TOOLS["bcftools"], "query", "-f", "%CHROM\t%POS\n", str(path)], text=True)
        self.assertNotIn("chr3\t", body)
        self.assertEqual(len(body.splitlines()), 6)
        self.assertEqual(original, (stale.read_bytes(), Path(str(stale) + ".tbi").read_bytes()))

    def test_unsupported_only_input_fails_without_publishing_stale_shards(self):
        folder, path, gnomad, _, _ = self.fixture("unsupported", stale=True, unsupported_only=True)
        before = path.read_bytes(), Path(str(path) + ".tbi").read_bytes()
        result = self.run_annotation(path, gnomad)
        self.assertNotEqual(result.returncode, 0)
        self.assertIn("not supported", result.stderr)
        self.assertEqual(before, (path.read_bytes(), Path(str(path) + ".tbi").read_bytes()))


FAILURES = {
    "query": 41, "split-view": 42, "split-sort": 43, "split-tabix": 44,
    "final-tabix": 45, "split-validation": 10, "annotated-validation": 47,
    "final-validation": 48, "merge-concat": 49, "merge-sort": 50,
    "parallel-dry": 51, "parallel-run": 52, "joblog-missing": 1,
    "joblog-empty": 1, "joblog-malformed": 2, "joblog-incomplete": 2,
    "joblog-nonzero": 17, "joblog-signal": 137,
}


def failure_test(failure, expected):
    def test(self):
        self.assert_failure(failure, expected)
    return test


for _failure, _expected in FAILURES.items():
    setattr(AnnotationFailures, "test_failure_" + _failure.replace("-", "_"),
            failure_test(_failure, _expected))


if __name__ == "__main__":
    unittest.main()
