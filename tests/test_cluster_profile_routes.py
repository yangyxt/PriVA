"""Current cluster profiles must not reuse another VEP release's readiness."""
import contextlib
import importlib.util
import io
import json
from pathlib import Path
import sys
import tempfile
import unittest
from unittest.mock import patch

import yaml

ROOT = Path(__file__).resolve().parents[1]


def load_script(name):
    spec = importlib.util.spec_from_file_location(name, ROOT / "scripts" / f"{name}.py")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


prepare = load_script("prepare_cluster_config")
readiness = load_script("check_cluster_readiness")


class RequestedClusterRelease(unittest.TestCase):
    def transform(self, config):
        return prepare.transform(config, Path("/source"), Path("/cache"), Path("/data"), 4)

    def test_requested_116_is_preserved_without_mutating_template(self):
        template = {"vep_installed_version": "116", "custom_setting": "retain"}
        result = self.transform(template)
        self.assertEqual(result["vep_installed_version"], "116")
        self.assertEqual(result["custom_setting"], "retain")
        self.assertEqual(template, {"vep_installed_version": "116", "custom_setting": "retain"})

    def test_integer_major_is_serialized_as_release_string(self):
        self.assertEqual(self.transform({"vep_installed_version": 116})["vep_installed_version"], "116")

    def test_explicit_historical_113_template_remains_113(self):
        self.assertEqual(self.transform({"vep_installed_version": "113"})["vep_installed_version"], "113")

    def test_missing_release_has_no_silent_default(self):
        with self.assertRaisesRegex(ValueError, "numeric VEP major release"):
            self.transform({})

    def test_invalid_major_is_rejected(self):
        for value in (None, "", "116.3", "release116", 116.0, True, "0"):
            with self.subTest(value=value), self.assertRaisesRegex(ValueError, "numeric VEP major release"):
                self.transform({"vep_installed_version": value})


class ReleaseAwareReadiness(unittest.TestCase):
    def setUp(self):
        temporary = tempfile.TemporaryDirectory(prefix="cluster_profile_", dir=ROOT / "tests")
        self.addCleanup(temporary.cleanup)
        self.folder = Path(temporary.name)

    def report(self, assembly="hg19", receipt_version="116.3", receipt_status="software_checked"):
        state = self.folder / "state"
        state.mkdir(exist_ok=True)
        receipt = state / "environment.ready.json"
        if receipt_version is not None:
            receipt.write_text(json.dumps({"status": receipt_status,
                                           "packages": {"ensembl-vep": receipt_version}}))
        elif receipt.exists():
            receipt.unlink()
        config = {
            "assembly": assembly, "vep_installed_version": "116",
            "vep_cache_dir": str(self.folder / "cache"),
            "gene_mechanism_builder_script": str(self.folder / "source/scripts/builder.py"),
            "clinvar_gene_stat": "clinvar.pkl", "interpro_mapping_pickle": "interpro.pkl",
            "loeuf_prescore": "loeuf.tsv", "repeat_region_file_name": "repeats.tsv",
        }
        config_path = self.folder / "config.yaml"
        config_path.write_text(yaml.safe_dump(config))
        output = self.folder / "report.json"
        command = ["check_cluster_readiness.py", "--config", str(config_path),
                   "--state-dir", str(state), "--output", str(output)]
        with patch.object(sys, "argv", command), contextlib.redirect_stdout(io.StringIO()):
            readiness.main()
        return json.loads(output.read_text())

    def test_113_receipt_cannot_qualify_116_configuration(self):
        report = self.report(receipt_version="113.2")
        self.assertFalse(report["software_qualified"])
        self.assertEqual(report["completed_stage_receipts"]["environment.ready.json"]
                         ["packages"]["ensembl-vep"], "113.2")

    def test_matching_116_receipt_qualifies_software_only(self):
        report = self.report()
        self.assertTrue(report["software_qualified"])
        self.assertEqual(report["status"], "incomplete_installation")
        self.assertIn("Shared PriVA/CADD official VEP116", [item["resource"] for item in report["blockers"]])
        self.assertIn("Family acceptance", [item["resource"] for item in report["blockers"]])

    def test_unchecked_matching_receipt_does_not_qualify_software(self):
        self.assertFalse(self.report(receipt_status="installing")["software_qualified"])

    def test_missing_receipt_does_not_qualify_software(self):
        self.assertFalse(self.report(receipt_version=None)["software_qualified"])

    def test_approved_hg19_mavedb_omission_is_not_a_download_blocker(self):
        report = self.report()
        self.assertEqual(report["approved_omissions"], ["MaveDB GRCh37 (explicit user approval)"])
        self.assertFalse(any("MaveDB" in item["resource"] for item in report["blockers"]))

    def test_grch38_reports_v2_and_retains_mavedb_requirement(self):
        report = self.report(assembly="GRCh38")
        labels = [item["resource"] for item in report["blockers"]]
        self.assertIn("SpliceAI2", labels)
        self.assertNotIn("SpliceAI v1", labels)
        self.assertEqual(report["approved_omissions"], [])
        resource = next(item for item in report["resources"] if item["key"] == "mavedb_file")
        self.assertEqual(resource["status"], "not_configured")

    def test_missing_or_malformed_package_identity_does_not_qualify(self):
        for receipt in ({}, {"status": "software_checked"},
                        {"status": "software_checked", "packages": []},
                        {"status": "software_checked", "packages": {"ensembl-vep": "116-unverified"}}):
            with self.subTest(receipt=receipt):
                self.assertFalse(readiness.environment_receipt_matches(receipt, "116"))

    def test_hg19_and_grch37_use_grch37_cache(self):
        for assembly in ("hg19", "GRCh37"):
            with self.subTest(assembly=assembly):
                report = self.report(assembly=assembly)
                entry = next(item for item in report["resources"] if item["key"] == "vep_cache_info")
                self.assertEqual(Path(entry["path"]).relative_to(self.folder).as_posix(),
                                 "cache/homo_sapiens_merged/116_GRCh37/info.txt")

    def test_hg38_and_grch38_use_grch38_cache(self):
        cache_info = self.folder / "cache/homo_sapiens_merged/116_GRCh38/info.txt"
        cache_info.parent.mkdir(parents=True)
        cache_info.write_text("fixture cache metadata\n")
        for assembly in ("hg38", "GRCh38"):
            with self.subTest(assembly=assembly):
                report = self.report(assembly=assembly)
                entry = next(item for item in report["resources"] if item["key"] == "vep_cache_info")
                self.assertEqual(entry["path"], str(cache_info))
                self.assertEqual(entry["status"], "present_not_scientifically_qualified")


if __name__ == "__main__":
    unittest.main()
