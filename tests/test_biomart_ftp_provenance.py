"""Distinct FTP116 provenance/reuse contracts; no scientific profile selection."""
from copy import deepcopy
import importlib.util
import json
import os
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

ROOT = Path(__file__).resolve().parents[1]


def load(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


B = load("ftp_test_biomart", ROOT / "scripts/biomart_nm_to_enst.py")
P = load("ftp_test_prerequisites", ROOT / "scripts/installer_prerequisites.py")
F = B.FTP


class FTPProvenance(unittest.TestCase):
    def setUp(self):
        folder = Path(os.environ.get("PRIVA_TEST_TMPDIR", tempfile.gettempdir())) / "priva-biomart-ftp-provenance-tests"
        folder.mkdir(parents=True, exist_ok=True)
        self.temp = tempfile.TemporaryDirectory(dir=folder)
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)
        self.path = self.root / F.BASENAME
        self.path.write_text("\t".join(B.HEADERS) + "\n"
                             "ENSG1\tENST2\tA\tP00001\tNM_1\n"
                             "ENSG1\tENST1\tA\tP00001\tNM_1\n"
                             "ENSG1\tENST2\tA\tP00001\tNM_1\n"
                             "ENSG2\tENST3\t\t\t\n")
        for key, value in (("EXPECTED_BYTES", self.path.stat().st_size), ("EXPECTED_SHA256", B.sha256(self.path)),
                           ("EXPECTED_STATISTICS", B.validate_table(self.path))):
            manager = patch.object(F, key, value)
            manager.start(); self.addCleanup(manager.stop)
        manager = patch.object(P, "BIOMART", B)
        manager.start(); self.addCleanup(manager.stop)
        self.metadata = B.provenance_path(self.path)

    def record(self):
        return F.make_provenance(self.path, B.validate_table, "registered_existing_exact_snapshot")

    def test_exact_existing_table_registers_and_reuses_offline_without_raw_export(self):
        original = self.path.read_bytes()
        with patch.object(F, "urlopen", side_effect=AssertionError("Unexpected network")):
            first = B.acquire(self.path, Path(str(self.path) + ".biomart.raw.tsv"), source_id=F.SOURCE_ID)
            before = self.metadata.read_bytes(), self.metadata.stat().st_mtime_ns
            second = B.acquire(self.path, source_id=F.SOURCE_ID)
        self.assertEqual(first, second)
        self.assertEqual(before, (self.metadata.read_bytes(), self.metadata.stat().st_mtime_ns))
        self.assertEqual(self.path.read_bytes(), original)
        self.assertFalse(Path(str(self.path) + ".biomart.raw.tsv").exists())
        self.assertFalse((self.root / ".biomart116_ftp").exists())
        self.assertFalse(first["online_query_executed"])
        self.assertNotIn("completion_stamp_verified", first)
        self.assertFalse(first["scientific_qualified"])
        self.assertFalse(first["production_selected"])
        self.assertIsNotNone(B.mapping_identity(self.path)["provenance_sha256"])

    def test_all_method_release_checksum_order_and_multiplicity_fields_are_verified(self):
        base = self.record()
        changes = ({"assembly": "GRCh37.p13"}, {"mart_release": 113},
                   {"acquisition_method": "online_biomart"}, {"source_database": "ensembl_mart_113"},
                   {"query_sql": F.QUERY + " DISTINCT"}, {"order_policy": "sorted/deduplicated"},
                   {"multiplicity_policy": "drop duplicates"}, {"scientific_qualified": True},
                   {"production_selected": 0}, {"completion_stamp_verified": True},
                   {"query_xml": "pretend service response"}, {"schema_sha256": "0" * 64})
        for fields in changes:
            with self.subTest(fields=fields):
                with self.assertRaises(ValueError):
                    F.verify_record(self.path, {**base, **fields}, B.validate_table, "hg38")
        wrong = deepcopy(base)
        wrong["raw_sources"][0]["sha256"] = "0" * 64
        with self.assertRaises(ValueError):
            F.verify_record(self.path, wrong, B.validate_table, "hg38")
        wrong = deepcopy(base)
        wrong["statistics"]["duplicate_nm_enst_rows"] = 0
        with self.assertRaisesRegex(ValueError, "statistics"):
            F.verify_record(self.path, wrong, B.validate_table, "hg38")
        with self.assertRaisesRegex(ValueError, "producer identity"):
            F.verify_record(self.path, {**base, "registration_producer_sha256": None}, B.validate_table, "hg38")

    def test_fresh_projection_operation_does_not_imply_new_raw_downloads(self):
        record = F.make_provenance(self.path, B.validate_table, "projected_from_verified_ftp_sources")
        self.assertEqual(F.verify_record(self.path, record, B.validate_table, "hg38")["operation"], "projected_from_verified_ftp_sources")
        with self.assertRaisesRegex(ValueError, "Unknown FTP resource operation"):
            F.verify_record(self.path, {**record, "operation": "downloaded_projected_and_verified"}, B.validate_table, "hg38")

    def test_reviewed_table_order_and_duplicates_are_pinned_even_with_same_pairs(self):
        record = self.record()
        lines = self.path.read_text().splitlines(keepends=True)
        self.path.write_text(lines[0] + lines[2] + lines[1] + "".join(lines[3:]))
        with self.assertRaisesRegex(ValueError, "SHA256"):
            F.verify_record(self.path, record, B.validate_table, "hg38")

    def test_wrong_existing_bytes_and_provenance_are_never_replaced(self):
        self.path.write_text("preserve arbitrary existing data")
        with patch.object(F, "urlopen", side_effect=AssertionError("Unexpected network")):
            with self.assertRaises(ValueError):
                B.acquire(self.path, source_id=F.SOURCE_ID)
        self.assertEqual(self.path.read_text(), "preserve arbitrary existing data")
        self.assertFalse(self.metadata.exists())

    def test_invalid_existing_sidecar_is_not_recovered_by_silent_re_registration(self):
        original = self.path.read_bytes()
        self.metadata.write_text('{"source_id":"unreviewed"}\n')
        before = self.metadata.read_bytes()
        with self.assertRaises(ValueError):
            B.acquire(self.path, source_id=F.SOURCE_ID)
        self.assertEqual(self.path.read_bytes(), original)
        self.assertEqual(self.metadata.read_bytes(), before)

    def test_raw_source_failure_never_publishes_mapping_or_provenance(self):
        self.path.unlink()
        with patch.object(F, "acquire", side_effect=ValueError("Publisher checksum changed")):
            with self.assertRaisesRegex(ValueError, "checksum changed"):
                B.acquire(self.path, source_id=F.SOURCE_ID)
        self.assertFalse(self.path.exists())
        self.assertFalse(self.metadata.exists())

    def test_orphan_provenance_blocks_fresh_fetch(self):
        self.path.unlink()
        self.metadata.write_text("{}")
        before = self.metadata.read_bytes()
        with self.assertRaisesRegex(ValueError, "existing provenance"):
            B.acquire(self.path, source_id=F.SOURCE_ID)
        self.assertEqual(self.metadata.read_bytes(), before)

    def test_37_and_38_defaults_are_unchanged_and_116_is_explicit(self):
        self.assertEqual(B.source_for_assembly("hg19")["source_id"], B.SOURCE_ID)
        self.assertEqual(B.source_for_assembly("hg38")["source_id"], B.GRCH38_SOURCE_ID)
        self.assertEqual(P.nm_to_enst_path({"assembly": "hg38", "nm_to_enst_tsv": str(self.path)}, ROOT), self.path)
        with self.assertRaisesRegex(ValueError, "not compatible"):
            P.nm_to_enst_path({"assembly": "hg19", "nm_to_enst_tsv": str(self.path)}, ROOT)
        with self.assertRaisesRegex(ValueError, "cannot be used"):
            F.verify_record(self.path, self.record(), B.validate_table, "hg19")

    def test_all_missing_explicit_inputs_require_separate_explicit_acquisition(self):
        self.path.unlink()
        config = {"assembly": "hg38", "nm_to_enst_tsv": str(self.path)}
        with self.assertRaisesRegex(ValueError, "restore that selected input"):
            P.mapping_fetch_source(config, ROOT)
        with self.assertRaisesRegex(ValueError, "not compatible"):
            P.mapping_fetch_source({**config, "assembly": "hg19"}, ROOT)
        with self.assertRaisesRegex(ValueError, "restore that selected input"):
            P.mapping_fetch_source({**config, "nm_to_enst_tsv": str(self.root / "custom.tsv")}, ROOT)

    def test_main_preflight_rejects_missing_explicit_ftp_and_remains_read_only(self):
        self.path.unlink()
        config = {"assembly": "hg38", "nm_to_enst_tsv": str(self.path)}
        with patch.object(P, "validate_clinvar_stats"):
            with self.assertRaisesRegex(ValueError, "Missing, empty, or unreadable"):
                P.validate_main(config, ROOT, allow_missing_bundled_mapping=True)
            with self.assertRaises(ValueError):
                P.validate_main({**config, "nm_to_enst_tsv": str(self.root / "unknown.tsv")}, ROOT, allow_missing_bundled_mapping=True)
        self.assertFalse(self.path.exists())
        self.assertFalse(self.metadata.exists())

    def test_logical_sidecar_and_assembly_binding_survive_annex_style_payload_links(self):
        payload = self.root / "objects" / "opaque-payload"
        payload.parent.mkdir()
        self.path.rename(payload)
        self.path.symlink_to(payload)
        B.acquire(self.path, source_id=F.SOURCE_ID)
        config = {"assembly": "hg38", "nm_to_enst_tsv": str(self.path)}
        with patch.dict(os.environ, {"PRIVA_INSTALL_INSPECT_ROOTS": str(self.root)}):
            self.assertEqual(P.validate_mapping(config, ROOT), self.path)
            self.assertIsNotNone(B.mapping_identity(self.path)["provenance_sha256"])
            with self.assertRaises(ValueError):
                P.validate_mapping({**config, "assembly": "hg19"}, ROOT)

    def test_ftp_source_never_uses_an_online_registry(self):
        with patch.object(B, "_request", side_effect=AssertionError("Unexpected online BioMart")):
            with self.assertRaisesRegex(ValueError, "not an online"):
                B.source_metadata(source_id=F.SOURCE_ID)


if __name__ == "__main__":
    unittest.main(verbosity=2)
