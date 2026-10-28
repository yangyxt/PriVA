"""Date/ownership/download and original-stage refresh contracts; no real network."""
from __future__ import annotations

from datetime import date, datetime, timezone
import gzip
import hashlib
import io
import json
import os
from pathlib import Path
import pickle
import subprocess
import sys
import tempfile
import unittest
from unittest import mock


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
import clinvar_gene_refresh as refresh


def vcf_bytes(stamp="20260430"):
    return (f"##fileformat=VCFv4.2\n##fileDate={stamp}\n"
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n").encode()


def write_vcf(path, stamp="20260430"):
    path.parent.mkdir(parents=True, exist_ok=True)
    with gzip.open(path, "wb") as handle:
        handle.write(vcf_bytes(stamp))


def write_stat(path):
    path.parent.mkdir(parents=True, exist_ok=True)
    with gzip.open(path, "wb") as handle:
        pickle.dump({"GENE": {"ensembl_id": "ENSG1", "clinvar_pathogenic": set(),
                              "large_aachange": set(), "small_aachange": set(),
                              "putative_nmd_variants": set()}}, handle)


class ClinVarGeneRefreshTests(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory(prefix="priva_clinvar_refresh_",
            dir=os.environ.get("PRIVA_TEST_TMPDIR", os.environ.get("TMPDIR")))
        self.addCleanup(self.temporary.cleanup)
        self.directory = Path(self.temporary.name)
        self.source = self.directory / "source"
        self.baseline = self.source / "data/ClinVar/clinvar_2star_stats.pkl.gz"
        self.baseline.parent.mkdir(parents=True)
        self.baseline.write_bytes(b"opaque nonempty bundled bytes; never unpickle or hash for reuse")
        (self.source / "cache_resources.json").write_text(json.dumps({"schema_version": 1, "resources": {
            "clinvar_gene_statistics": {"refresh": {"threshold_months": 3, "automatic": True,
                "policy": "latest_filedate_vs_bundled_mtime"}},
            "clinvar_vcf": {"source": {"latest_url_template": refresh.NCBI + "/pub/clinvar/vcf_{assembly}/clinvar.vcf.gz"}},
            "clinvar_variant_summary": {"source": {"url": refresh.SUMMARY_URL}},
            "mane_select": {"source": {"url": refresh.MANE_URL}},
        }}))
        self.baseline_time = datetime(2026, 1, 31, 20, 15, tzinfo=timezone.utc).timestamp()
        os.utime(self.baseline, (self.baseline_time, self.baseline_time))
        for name in ("clinvar_stat_variants.py", "combine_annotations.py", "annotation_vcf.sh"):
            path = self.source / "scripts" / name
            path.parent.mkdir(parents=True, exist_ok=True)
            path.write_text("frozen fixture scientific source\n")
        self.cache = self.directory / "cache"
        self.config = self.directory / "config.json"
        self.settings = {"assembly": "hg38", "vep_installed_version": "116",
                         "clinvar_vcf_dir": str(self.cache), "threads": 1,
                         "clinvar_gene_stat": "clinvar_2star_stats.pkl.gz"}
        self.save_config()
        self.latest = {"url": "https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/clinvar.vcf.gz",
                       "file_date": "2026-04-30", "observed_utc": "2026-05-01T00:00:00+00:00",
                       "publisher_md5": "1" * 32}
        self.metadata = mock.patch.object(refresh, "latest_metadata", return_value=self.latest).start()
        self.addCleanup(mock.patch.stopall)
        self.calls = []

    def save_config(self):
        self.config.write_text(json.dumps(self.settings))

    def plan(self):
        return refresh.plan_refresh(self.config, self.source)

    def acquire_fixture(self, plan):
        folder = Path(plan["build_dir"]) / "inputs"
        folder.mkdir(parents=True, exist_ok=True)
        vcf = folder / "dated.vcf.gz"
        write_vcf(vcf, plan["latest"]["file_date"].replace("-", ""))
        Path(str(vcf) + ".tbi").write_bytes(b"source index")
        sources = {"vcf": refresh.identity(vcf)}
        for name in ("summary", "mane"):
            path = folder / (name + ".txt.gz")
            path.write_bytes(gzip.compress(("pinned " + name).encode()))
            sources[name] = refresh.identity(path)
        return {"sources": sources, "mane_release": "1.4"}

    def original_fixture(self, source_root, function, *arguments):
        self.assertEqual(source_root, self.source)
        self.calls.append((function, arguments))
        current = refresh.load_config(Path(arguments[0]))
        if function == "ClinVar_VCF_deploy":
            path = refresh.selected_vcf(current, Path(current["clinvar_vcf_dir"]), "hg38")
            write_vcf(path)
            Path(str(path) + ".tbi").write_bytes(b"completed annotation index")
        elif function == "ClinVar_AA_stat":
            for key in ("clinvar_aa_stat", "clinvar_splice_stat"):
                path = Path(current[key])
                path.parent.mkdir(parents=True, exist_ok=True)
                path.write_bytes(b"original AA/splice output")
        elif function == "ClinVar_Gene_stat_regenerate":
            write_stat(Path(arguments[3]))

    def install(self, original=None):
        with mock.patch.object(refresh, "validate_main") as preflight, \
             mock.patch.object(refresh, "acquire_inputs", side_effect=self.acquire_fixture) as acquired, \
             mock.patch.object(refresh, "run_original", side_effect=original or self.original_fixture):
            value = refresh.install_refresh(self.config, self.source)
        return value, preflight, acquired

    def test_calendar_months_include_boundary_and_clamp_month_ends(self):
        for start, expected in ((date(2026, 1, 31), date(2026, 4, 30)),
                                (date(2023, 11, 30), date(2024, 2, 29)),
                                (date(2026, 11, 30), date(2027, 2, 28)),
                                (date(2026, 10, 17), date(2027, 1, 17))):
            with self.subTest(start=start):
                self.assertEqual(refresh.add_calendar_months(start), expected)
        self.assertEqual(self.plan()["action"], "rebuild")
        self.latest["file_date"] = "2026-04-29"
        self.assertEqual(self.plan()["action"], "preserve")

    def test_bundle_is_metadata_only_and_file_date_controls_age(self):
        with mock.patch.object(refresh, "digest", wraps=refresh.digest) as hashed:
            plan = self.plan()
        self.assertNotIn(self.baseline, [call.args[0] for call in hashed.call_args_list])
        self.assertEqual(plan["baseline"]["mtime_utc"], "2026-01-31T20:15:00+00:00")
        selected = self.directory / "explicit.vcf.gz"
        write_vcf(selected, "20260429")
        future = datetime(2030, 1, 1, tzinfo=timezone.utc).timestamp()
        os.utime(selected, (future, future))
        self.settings["clinvar_vcf"] = str(selected)
        self.save_config()
        with self.assertRaisesRegex(ValueError, "Explicit ClinVar VCF.*2026-04-29.*2026-04-30"):
            self.plan()

    def test_both_assemblies_probe_own_official_release(self):
        for assembly, expected in (("hg19", "GRCh37"), ("GRCh37", "GRCh37"),
                                   ("hg38", "GRCh38"), ("GRCh38", "GRCh38")):
            with self.subTest(assembly=assembly):
                self.settings["assembly"] = assembly
                self.save_config()
                self.assertEqual(self.plan()["assembly"], expected)
                self.metadata.assert_called_with(expected, url_template=refresh.NCBI + "/pub/clinvar/vcf_{assembly}/clinvar.vcf.gz")

    def test_old_managed_vcf_does_not_hide_latest_publication(self):
        old = self.cache / "clinvar.hg38.vep.vcf.gz"
        write_vcf(old, "20260101")
        self.settings["clinvar_vcf"] = str(old)
        self.save_config()
        before = old.read_bytes()
        plan = self.plan()
        self.assertFalse(plan["explicit_vcf"])
        self.assertIn("GRCh38/20260430/vep116", plan["output_path"])
        self.assertEqual(old.read_bytes(), before)

    def test_explicit_newer_or_undated_vcf_is_not_silently_replaced(self):
        explicit = self.directory / "chosen.vcf.gz"
        self.settings["clinvar_vcf"] = str(explicit)
        self.save_config()
        for payload in (vcf_bytes("20260501"), b"##fileformat=VCFv4.2\n#CHROM\tPOS\n"):
            compressed = gzip.compress(payload)
            explicit.write_bytes(compressed)
            with self.assertRaises(ValueError):
                self.plan()
            self.assertEqual(explicit.read_bytes(), compressed)

    def test_explicit_selected_stat_path_is_preserved_without_a_freshness_claim(self):
        explicit = self.directory / "selected directory" / "chosen.pkl.gz"
        write_stat(explicit)
        self.settings["clinvar_gene_stat"] = str(explicit)
        self.save_config()
        before = self.config.read_bytes()
        plan = self.plan()
        self.assertEqual(plan["action"], "preserve")
        self.metadata.assert_not_called()
        self.assertEqual(self.config.read_bytes(), before)
        self.baseline.unlink()
        self.assertEqual(self.plan()["action"], "preserve")

    def test_success_uses_original_stage_order_and_receipt_prevents_repeat_rebuild(self):
        before = self.baseline.read_bytes()
        result, preflight, acquired = self.install()
        self.assertEqual([call[0] for call in self.calls],
                         ["ClinVar_VCF_deploy", "ClinVar_AA_stat", "ClinVar_Gene_stat_regenerate"])
        self.assertEqual(self.calls[0][1][1], "--source-vcf")
        preflight.assert_called_once()
        acquired.assert_called_once()
        selected = Path(refresh.load_config(self.config)["clinvar_gene_stat"])
        self.assertEqual(selected, Path(result["output_path"]))
        receipt = json.loads(refresh.receipt_path(selected).read_text())
        self.assertEqual(receipt["data_date"], "2026-04-30")
        self.assertEqual(receipt["status"], "complete")
        self.assertEqual(self.baseline.read_bytes(), before)
        self.assertEqual(self.baseline.stat().st_mtime, self.baseline_time)
        self.calls.clear()
        again, _, acquired = self.install()
        self.assertEqual(again["action"], "reuse")
        acquired.assert_not_called()
        self.assertEqual(self.calls, [])

    def test_missing_baseline_is_reconstructed_in_managed_cache_without_recreating_bundle(self):
        self.baseline.unlink()
        self.latest["file_date"] = "2026-01-01"  # No age boundary exists without a baseline.
        def original(source_root, function, *arguments):
            self.original_fixture(source_root, function, *arguments)
            if function == "ClinVar_VCF_deploy":
                current = refresh.load_config(Path(arguments[0]))
                write_vcf(Path(current["clinvar_vcf"]), "20260101")
        result, preflight, _ = self.install(original)
        self.assertTrue(result["baseline_missing"])
        self.assertFalse(self.baseline.exists())
        self.assertTrue(Path(result["output_path"]).is_relative_to(self.cache))
        preflight.assert_called_once_with(mock.ANY, self.source, allow_missing_bundled_clinvar=True)
        self.calls.clear()
        again, _, acquired = self.install()
        self.assertEqual(again["action"], "reuse")
        acquired.assert_not_called()
        self.assertEqual(self.calls, [])

    def test_empty_baseline_is_metadata_only_and_not_overwritten(self):
        self.baseline.write_bytes(b"")
        result, _, _ = self.install()
        self.assertTrue(result["baseline_missing"])
        self.assertEqual(self.baseline.read_bytes(), b"")
        receipt = json.loads(refresh.receipt_path(Path(result["output_path"])).read_text())
        self.assertIs(receipt["scientific_qualified"], False)

    def test_missing_explicit_statistics_remain_an_error(self):
        self.settings["clinvar_gene_stat"] = str(self.cache / "explicit-missing.pkl.gz")
        self.save_config()
        with self.assertRaises(OSError):
            self.plan()
        self.metadata.assert_not_called()

    def test_failed_missing_baseline_reconstruction_preserves_config_and_absence(self):
        self.baseline.unlink()
        before = self.config.read_bytes()
        def fail(source_root, function, *arguments):
            if function == "ClinVar_Gene_stat_regenerate":
                raise subprocess.CalledProcessError(17, function)
            self.original_fixture(source_root, function, *arguments)
        with self.assertRaises(subprocess.CalledProcessError):
            self.install(fail)
        self.assertEqual(self.config.read_bytes(), before)
        self.assertFalse(self.baseline.exists())
        self.assertFalse(refresh.receipt_path(Path(self.plan()["output_path"])).exists())

    def test_new_baseline_during_missing_baseline_build_prevents_selection(self):
        self.baseline.unlink()
        before = self.config.read_bytes()
        def original(source_root, function, *arguments):
            self.original_fixture(source_root, function, *arguments)
            if function == "ClinVar_Gene_stat_regenerate":
                self.baseline.write_bytes(b"concurrent bundled snapshot")
        with self.assertRaisesRegex(ValueError, "Bundled ClinVar metadata changed"):
            self.install(original)
        self.assertEqual(self.config.read_bytes(), before)

    def test_failed_generation_preserves_original_config_and_bundled_file(self):
        before_config = self.config.read_bytes()
        before_bundle = self.baseline.read_bytes()
        def fail(source_root, function, *arguments):
            if function == "ClinVar_Gene_stat_regenerate":
                raise subprocess.CalledProcessError(17, function)
            self.original_fixture(source_root, function, *arguments)
        with self.assertRaises(subprocess.CalledProcessError):
            self.install(fail)
        self.assertEqual(self.config.read_bytes(), before_config)
        self.assertEqual(self.baseline.read_bytes(), before_bundle)
        self.assertFalse(refresh.receipt_path(Path(self.plan()["output_path"])).exists())

    def test_reset_to_bundled_stat_with_identical_completed_vcf_reselects_ready_result(self):
        result, _, _ = self.install()
        self.settings = refresh.load_config(self.config)
        self.settings["clinvar_gene_stat"] = "clinvar_2star_stats.pkl.gz"
        self.save_config()
        self.calls.clear()
        reused, _, acquired = self.install()
        self.assertEqual(reused["action"], "reuse")
        acquired.assert_not_called()
        self.assertEqual(self.calls, [])
        self.assertEqual(refresh.load_config(self.config)["clinvar_gene_stat"], result["output_path"])

    def test_changed_profile_uses_new_workspace_and_preserves_failed_prior_annotation(self):
        def fail(source_root, function, *arguments):
            if function == "ClinVar_Gene_stat_regenerate":
                raise subprocess.CalledProcessError(17, function)
            self.original_fixture(source_root, function, *arguments)
        with self.assertRaises(subprocess.CalledProcessError):
            self.install(fail)
        folder = Path(self.plan()["build_dir"])
        plan_bytes = (folder / "refresh.plan.json").read_bytes()
        config_bytes = (folder / "install.config.yaml").read_bytes()
        annotation = folder / "clinvar.hg38.vep.vcf.gz"
        annotation_bytes = annotation.read_bytes()
        (self.source / "scripts/clinvar_stat_variants.py").write_text("different profile B\n")
        rebuilt, _, acquired = self.install()
        acquired.assert_called_once()
        self.assertNotEqual(Path(rebuilt["build_dir"]), folder)
        self.assertEqual((folder / "refresh.plan.json").read_bytes(), plan_bytes)
        self.assertEqual((folder / "install.config.yaml").read_bytes(), config_bytes)
        self.assertEqual(annotation.read_bytes(), annotation_bytes)

    def test_unknown_preexisting_annotation_workspace_is_not_blessed(self):
        plan = self.plan()
        folder = Path(plan["build_dir"])
        write_vcf(folder / "clinvar.hg38.vep.vcf.gz")
        with mock.patch.object(refresh, "validate_main"), \
             mock.patch.object(refresh, "acquire_inputs") as acquired:
            with self.assertRaisesRegex(ValueError, "no owned scientific-profile plan"):
                refresh.install_refresh(self.config, self.source)
        acquired.assert_not_called()
        self.assertFalse((folder / "refresh.plan.json").exists())

    def test_preflight_failure_precedes_network_and_any_acquisition(self):
        with mock.patch.object(refresh, "validate_main", side_effect=ValueError("missing selected input")), \
             mock.patch.object(refresh, "acquire_inputs") as acquired:
            with self.assertRaisesRegex(ValueError, "missing selected"):
                refresh.install_refresh(self.config, self.source)
        self.metadata.assert_not_called()
        acquired.assert_not_called()
        self.assertFalse(self.cache.exists())

    def test_network_failure_is_reported_without_cached_freshness_fallback(self):
        self.metadata.side_effect = ValueError("cannot determine latest official ClinVar")
        with self.assertRaisesRegex(ValueError, "latest official"):
            self.plan()
        self.assertFalse(self.cache.exists())

    def test_missing_or_malformed_resource_registry_fails_before_network(self):
        path = self.source / "cache_resources.json"
        for content in ("{}", '{"schema_version":1,"resources":{}}', "invalid JSON"):
            with self.subTest(content=content):
                path.write_text(content)
                with self.assertRaisesRegex(ValueError, "Cannot load ClinVar refresh policy"):
                    self.plan()
        path.unlink()
        with self.assertRaisesRegex(ValueError, "Cannot load ClinVar refresh policy"):
            self.plan()
        self.metadata.assert_not_called()

    def test_registry_can_explicitly_disable_automatic_refresh_without_a_network_claim(self):
        path = self.source / "cache_resources.json"
        policy = json.loads(path.read_text())
        policy["resources"]["clinvar_gene_statistics"]["refresh"]["automatic"] = False
        path.write_text(json.dumps(policy))
        self.assertEqual(self.plan()["action"], "preserve")
        self.metadata.assert_not_called()

    def test_changed_owned_stat_is_not_replaced_or_reused_silently(self):
        self.install()
        selected = Path(refresh.load_config(self.config)["clinvar_gene_stat"])
        selected.write_bytes(b"externally changed")
        with self.assertRaisesRegex(ValueError, "differ from their completed receipt"):
            self.plan()

    def test_changed_scientific_source_builds_new_same_date_cache_and_preserves_prior_bytes(self):
        previous, _, _ = self.install()
        old = Path(previous["output_path"])
        old_bytes = old.read_bytes()
        old_receipt = refresh.receipt_path(old).read_bytes()
        before = self.config.read_bytes()
        (self.source / "scripts/clinvar_stat_variants.py").write_text("different scientific operation\n")
        planned = self.plan()
        self.assertEqual(planned["action"], "rebuild")
        self.assertNotEqual(planned["output_path"], previous["output_path"])
        self.assertEqual(self.config.read_bytes(), before)
        rebuilt, _, _ = self.install()
        self.assertEqual(rebuilt["output_path"], planned["output_path"])
        self.assertEqual(old.read_bytes(), old_bytes)
        self.assertEqual(refresh.receipt_path(old).read_bytes(), old_receipt)
        self.assertEqual(self.plan()["action"], "reuse")

    def test_changed_selected_reference_identity_gets_distinct_same_date_destination(self):
        reference = self.directory / "reference.fa"
        reference.write_text(">chr1\nACGT\n")
        self.settings["ref_genome"] = str(reference)
        self.save_config()
        previous, _, _ = self.install()
        previous_mtime = reference.stat().st_mtime_ns
        reference.write_text(">chr1\nAGGT\n")
        os.utime(reference, ns=(previous_mtime, previous_mtime))
        planned = self.plan()
        self.assertEqual(planned["action"], "rebuild")
        self.assertNotEqual(planned["output_path"], previous["output_path"])

    def test_identical_hashed_source_with_changed_mtime_reuses_same_profile(self):
        previous, _, _ = self.install()
        path = self.source / "scripts/clinvar_stat_variants.py"
        timestamp = path.stat().st_mtime_ns + 1_000_000_000
        os.utime(path, ns=(timestamp, timestamp))
        planned = self.plan()
        self.assertEqual(planned["action"], "reuse")
        self.assertEqual(planned["output_path"], previous["output_path"])
        self.assertEqual(planned["scientific_profile_sha256"], previous["scientific_profile_sha256"])

    def test_same_profile_legacy_unsuffixed_completed_cache_remains_reusable(self):
        legacy = self.plan()
        folder = Path(legacy["build_dir"]).with_name("vep116")
        legacy.update(build_dir=str(folder), output_path=str(folder / "clinvar_2star_stats.pkl.gz"))
        with mock.patch.object(refresh, "plan_refresh", return_value=legacy):
            self.install()
        self.calls.clear()
        reused, _, acquired = self.install()
        self.assertEqual(reused["action"], "reuse")
        self.assertEqual(reused["build_dir"], str(folder))
        acquired.assert_not_called()
        self.assertEqual(self.calls, [])

    def test_unrelated_new_file_in_shared_cache_container_does_not_change_profile(self):
        container = self.directory / "shared plugin caches"
        container.mkdir()
        self.settings["vep_plugins_cachedir"] = str(container)
        self.save_config()
        self.install()
        (container / "unrelated later resource").write_text("not a selected scientific input")
        self.assertEqual(self.plan()["action"], "reuse")

    def test_missing_or_changed_annotation_aa_splice_artifacts_block_reuse(self):
        self.install()
        selected = Path(refresh.load_config(self.config)["clinvar_gene_stat"])
        receipt = json.loads(refresh.receipt_path(selected).read_text())
        for name in ("vcf", "vcf_index", "aa", "splice"):
            path = Path(receipt["completed_artifacts"][name]["path"])
            original = path.read_bytes()
            for changed in (None, b"different completed bytes"):
                with self.subTest(name=name, missing=changed is None):
                    if changed is None:
                        path.unlink()
                    else:
                        path.write_bytes(changed)
                    with self.assertRaisesRegex(ValueError, "Completed ClinVar artifact"):
                        self.plan()
                    path.write_bytes(original)

    def test_same_date_explicit_different_vcf_cannot_reselect_old_annotation(self):
        self.install()
        chosen = self.directory / "explicit same date.vcf.gz"
        write_vcf(chosen)
        self.settings = refresh.load_config(self.config)
        self.settings["clinvar_vcf"] = str(chosen)
        self.save_config()
        before = self.config.read_bytes()
        with self.assertRaisesRegex(ValueError, "Explicit ClinVar VCF identity differs"):
            self.plan()
        self.assertEqual(self.config.read_bytes(), before)

    def test_same_date_new_publisher_bytes_are_not_treated_as_the_completed_release(self):
        self.install()
        self.latest["publisher_md5"] = "2" * 32
        with self.assertRaisesRegex(ValueError, "Official ClinVar bytes changed"):
            self.plan()

    def test_full_main_uses_the_named_clinvar_refresh_entrypoint_before_alpha_derivatives(self):
        reference = self.directory / "reference.fa"
        reference.write_text(">chr1\nA\n")
        environment = dict(os.environ, TEST_CONFIG=str(self.config), TEST_FOLDER=str(self.directory),
                           TEST_REFERENCE=str(reference), PRIVA_SMOKE_SKIP_VEP="1")
        script = r'''
source "$1" >/dev/null || exit
log() { :; }
priva_install_basic_paths() { PRIVA_INSTALL_TMPDIR="$TEST_FOLDER"; }
conda_install_vep() { PRIVA_INSTALL_ENV_NAME=priva_acmg; }
priva_prepare_install_config() { :; }
read_yaml() { case "$2" in ref_genome) printf '%s' "$TEST_REFERENCE";; assembly) printf hg38;; *) printf '%s' "$TEST_FOLDER";; esac; }
reference_genome_install() { :; }
gnomAD_install() { :; }
AlphaMissense_anno_gnomAD() { :; }
InterPro_parsing() { :; }
ClinVar_install() { [[ $# == 2 && $1 == --config && $2 == "$TEST_CONFIG" ]] || return 91; printf 'CLINVAR_AUTO\n'; }
AlphaMissense_stat() { printf 'ALPHA_STAT\n'; }
AlphaMissense_pick_intolerant_domains() { :; }
prepare_pm1_regions() { :; }
ClinVar_patho_AF_stat() { printf 'CLINVAR_AF\n'; }
CADD_install() { :; }
CADD_runtime_install() { :; }
CADD_upstream_check() { :; }
ClinGen_deploy() { :; }
CDS_FASTA_install() { :; }
pext_install() { :; }
mechanism_resource_install() { :; }
priva_export_installed_environment() { :; }
main_install "$TEST_CONFIG"
'''
        result = subprocess.run(["bash", "-c", script, "fixture", str(ROOT / "scripts/install_utils.sh")],
                                env=environment, capture_output=True, text=True, timeout=20)
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertEqual(result.stdout.splitlines(), ["CLINVAR_AUTO", "ALPHA_STAT", "CLINVAR_AF"])


class ClinVarMetadataAndAcquisitionTests(unittest.TestCase):
    def test_remote_header_is_bounded_and_uses_file_date(self):
        for assembly in ("GRCh37", "GRCh38"):
            for stamp in ("20260430", "2026-04-30"):
                with self.subTest(assembly=assembly, stamp=stamp):
                    payload = gzip.compress(vcf_bytes(stamp))
                    def response(url, **kwargs):
                        return io.BytesIO((hashlib.md5(payload).hexdigest() + " source\n").encode() if url.endswith(".md5") else payload)
                    with mock.patch.object(refresh, "open_url", side_effect=response) as opened:
                        self.assertEqual(refresh.latest_metadata(assembly, url_template=refresh.NCBI + "/pub/clinvar/vcf_{assembly}/clinvar.vcf.gz")["file_date"], "2026-04-30")
                    self.assertIn(mock.call(f"{refresh.NCBI}/pub/clinvar/vcf_{assembly}/clinvar.vcf.gz", header_only=True), opened.call_args_list)

    def test_invalid_truncated_and_oversized_headers_fail_closed(self):
        payloads = [b"<html>unavailable</html>", gzip.compress(vcf_bytes())[:15],
                    gzip.compress(b"##fileformat=VCFv4.2\n##large=" + b"x" * 65537 + b"\n" + vcf_bytes())]
        for payload in payloads:
            with self.subTest(size=len(payload)):
                with mock.patch.object(refresh, "publisher_md5", return_value="0" * 32), \
                     mock.patch.object(refresh, "open_url", return_value=io.BytesIO(payload)):
                    with self.assertRaisesRegex(ValueError, "Cannot determine latest official ClinVar"):
                        refresh.latest_metadata("GRCh38", url_template=refresh.NCBI + "/pub/clinvar/vcf_{assembly}/clinvar.vcf.gz")

    def test_acquisition_pins_summary_independently_and_reuses_completed_bytes(self):
        with tempfile.TemporaryDirectory(dir=os.environ.get("PRIVA_TEST_TMPDIR", os.environ.get("TMPDIR"))) as temporary:
            folder = Path(temporary)
            vcf = gzip.compress(vcf_bytes())
            summary = gzip.compress(b"#AlleleID\tAssembly\tName\n")
            mane = gzip.compress(b"RefSeq_nuc\tEnsembl_nuc\n")
            vcf_url = refresh.NCBI + "/pub/clinvar/vcf_GRCh38/clinvar.vcf.gz"
            payloads = {vcf_url: vcf, vcf_url + ".tbi": b"index", refresh.SUMMARY_URL: summary, refresh.MANE_URL: mane}
            fetched = []
            def network(url, **kwargs):
                fetched.append((url, kwargs.get("header_only", False)))
                if url.endswith(".md5"):
                    return io.BytesIO((hashlib.md5(payloads[url[:-4]]).hexdigest() + "  source\n").encode())
                return io.BytesIO(payloads[url])
            plan = {"build_dir": str(folder), "assembly": "GRCh38", "explicit_vcf": False,
                    "latest": {"file_date": "2026-04-30", "publisher_md5": hashlib.md5(vcf).hexdigest()}, "resource_policy": {
                        "latest_url_template": refresh.NCBI + "/pub/clinvar/vcf_{assembly}/clinvar.vcf.gz",
                        "summary_url": refresh.SUMMARY_URL, "mane_url": refresh.MANE_URL}}
            with mock.patch.object(refresh, "open_url", side_effect=network):
                first = refresh.acquire_inputs(plan)
                fetched.clear()
                again = refresh.acquire_inputs(plan)
            self.assertEqual(first["sources"], again["sources"])
            self.assertTrue(all(url.endswith(".md5") or header_only for url, header_only in fetched))
            self.assertIn("independent publisher MD5", first["summary_release_identity"])
            self.assertEqual(first["mane_release"], "1.4")
            self.assertTrue((folder / "inputs.ready.json").exists())

    def test_download_checksum_failure_does_not_publish_or_overwrite_existing_file(self):
        with tempfile.TemporaryDirectory(dir=os.environ.get("PRIVA_TEST_TMPDIR", os.environ.get("TMPDIR"))) as temporary:
            path = Path(temporary) / "source.gz"
            path.write_bytes(b"retained earlier file")
            with mock.patch.object(refresh, "open_url", return_value=io.BytesIO(gzip.compress(b"new wrong bytes"))):
                with self.assertRaisesRegex(ValueError, "differ from publisher MD5"):
                    refresh.acquire_file("https://ftp.ncbi.nlm.nih.gov/source.gz", path, expected_md5="0" * 32)
            self.assertEqual(path.read_bytes(), b"retained earlier file")
            self.assertFalse(Path(str(path) + ".source.json").exists())


if __name__ == "__main__":
    unittest.main()
