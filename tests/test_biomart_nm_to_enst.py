"""Pinned BioMart export, immutable reuse and derived-output provenance contracts."""
import contextlib
import fcntl
import hashlib
import importlib.util
import io
import json
import os
from pathlib import Path
import tempfile
import unittest
from unittest import mock


ROOT = Path(__file__).resolve().parents[1]
SPEC = importlib.util.spec_from_file_location("biomart_mapping", ROOT / "scripts/biomart_nm_to_enst.py")
MAPPING = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(MAPPING)
REGISTRY = ('<MartRegistry><MartURLLocation name="ENSEMBL_MART_ENSEMBL" '
            'database="ensembl_mart_116" displayName="Ensembl Genes 116" '
            'host="grch37.ensembl.org" path="/biomart/martservice" '
            'serverVirtualSchema="default" /></MartRegistry>')
DATASETS = "TableSet\thsapiens_gene_ensembl\tHuman genes (GRCh37.p13)\t1\tGRCh37.p13\t200\t50000\tdefault\t2022-01-29 13:33:38\n"
BODY = ("\t".join(MAPPING.HEADERS) + "\n"
        "ENSG000001\tENST000001\tGENE\tP00001\tNM_000001\n"
        "ENSG000001\tENST000001\tGENE\tP00001\tNM_000001\n"
        "ENSG000001\tENST000002\tGENE\tP00002\tNM_000001\n"
        "ENSG000002\tENST000003\t\t\tNC_000016\n"
        "ENSG000003\tENST000004\t\t\t\n").encode()
REGISTRY38 = REGISTRY.replace("116", "113").replace("grch37.ensembl.org", "oct2024.archive.ensembl.org")
DATASETS38 = DATASETS.replace("GRCh37.p13", "GRCh38.p14")


class BioMartMappingTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory(prefix="biomart_mapping_", dir=os.environ.get("PRIVA_TEST_TMPDIR", os.environ.get("TMPDIR")))
        self.addCleanup(self.temp.cleanup)
        self.folder = Path(self.temp.name)
        self.output = self.folder / MAPPING.BASENAME
        self.raw = self.folder / "source.raw.tsv"
        self.metadata = MAPPING.parse_metadata(REGISTRY, DATASETS)
        snapshot = mock.patch.multiple(MAPPING, EXPECTED_SNAPSHOT_BYTES=len(BODY),
                                       EXPECTED_SNAPSHOT_SHA256=hashlib.sha256(BODY).hexdigest())
        snapshot.start()
        self.addCleanup(snapshot.stop)

    def acquire(self, body=BODY + b"[success]\n", metadata=None):
        source = metadata if metadata is not None else self.metadata
        with mock.patch.object(MAPPING, "source_metadata", return_value=source), \
                mock.patch.object(MAPPING, "_request", return_value=io.BytesIO(body)), \
                contextlib.redirect_stdout(io.StringIO()):
            return MAPPING.acquire(self.output, self.raw, attempts=1)

    def test_complete_export_preserves_exact_order_duplicates_and_ignored_non_nm(self):
        record = self.acquire()
        self.assertEqual(self.output.read_bytes(), BODY)
        self.assertEqual(self.raw.read_bytes(), BODY + b"[success]\n")
        self.assertEqual(record["statistics"]["nm_rows"], 3)
        self.assertEqual(record["statistics"]["duplicate_full_rows"], 1)
        self.assertEqual(record["statistics"]["duplicate_nm_enst_rows"], 1)
        self.assertEqual(record["statistics"]["nm_with_multiple_enst"], 1)
        self.assertEqual(record["statistics"]["non_nm_refseq_accessions"], ["NC_000016"])
        self.assertTrue(record["completion_stamp_verified"])
        self.assertTrue(Path(str(self.raw) + ".metadata-before.json").is_file())
        self.assertEqual(MAPPING.verify_pinned(self.output, "hg19")["sha256"], MAPPING.sha256(self.output))

    def test_completed_output_reuses_offline_and_rejects_changed_bytes(self):
        self.acquire()
        with mock.patch.object(MAPPING, "source_metadata", side_effect=AssertionError("network")), \
                contextlib.redirect_stdout(io.StringIO()):
            MAPPING.acquire(self.output, self.raw, attempts=1)
        self.output.write_bytes(BODY.replace(b"P00001", b"P99999"))
        with self.assertRaisesRegex(ValueError, "SHA256 mismatch"):
            MAPPING.verify_pinned(self.output)

    def test_partial_or_error_response_is_never_installed(self):
        for body in (BODY, b"Query ERROR: unavailable\n", BODY + b"[success]\nTRAILING ERROR\n"):
            with self.subTest(body=body[:20]):
                with self.assertRaises(ValueError):
                    self.acquire(body)
                self.assertFalse(self.output.exists())
                self.assertFalse(MAPPING.provenance_path(self.output).exists())

    def test_versioned_nm_is_rejected_without_rewriting(self):
        body = BODY.replace(b"NM_000001", b"NM_000001.4") + b"[success]\n"
        with self.assertRaisesRegex(ValueError, "bare IDs"):
            self.acquire(body)
        self.assertFalse(self.output.exists())
        self.assertEqual(Path(str(self.raw) + ".part").read_bytes(), body)

    def test_schema_failure_does_not_repeat_whole_export(self):
        with mock.patch.object(MAPPING, "source_metadata", return_value=self.metadata), \
                mock.patch.object(MAPPING, "_request", return_value=io.BytesIO(b"HTML error\n")) as request, \
                contextlib.redirect_stdout(io.StringIO()):
            with self.assertRaisesRegex(ValueError, "schema"):
                MAPPING.acquire(self.output, self.raw, attempts=3)
        self.assertEqual(request.call_count, 1)

    def test_mart_rollover_or_assembly_mismatch_is_explicit(self):
        with self.assertRaisesRegex(ValueError, "Pinned BioMart source changed"):
            MAPPING.parse_metadata(REGISTRY.replace("116", "117"), DATASETS)
        with self.assertRaisesRegex(ValueError, "dataset/assembly changed"):
            MAPPING.parse_metadata(REGISTRY, DATASETS.replace("GRCh37.p13", "GRCh38.p14"))
        self.acquire()
        with self.assertRaisesRegex(ValueError, "cannot be selected implicitly"):
            MAPPING.verify_pinned(self.output, "hg38")

    def test_metadata_change_during_export_prevents_publication(self):
        after = MAPPING.parse_metadata(REGISTRY, DATASETS.replace("2022-01-29", "2022-01-30"))
        with mock.patch.object(MAPPING, "source_metadata", side_effect=[self.metadata, after]), \
                mock.patch.object(MAPPING, "_request", return_value=io.BytesIO(BODY + b"[success]\n")), \
                contextlib.redirect_stdout(io.StringIO()):
            with self.assertRaisesRegex(ValueError, "changed during export"):
                MAPPING.acquire(self.output, self.raw, attempts=1)
        self.assertFalse(self.output.exists())

    def test_same_release_changed_bytes_or_order_cannot_publish_a_new_snapshot(self):
        lines = BODY.splitlines(keepends=True)
        reordered = b"".join([lines[0], lines[3], lines[1], lines[2], *lines[4:]])
        for changed in (BODY.replace(b"P00001", b"P99999"), reordered,
                        BODY + b"ENSG000009\tENST000009\t\t\tNM_000009\n"):
            with self.subTest(changed_sha256=hashlib.sha256(changed).hexdigest()):
                with self.assertRaisesRegex(ValueError, "default snapshot size/SHA256 mismatch"):
                    self.acquire(changed + b"[success]\n")
                self.assertFalse(self.output.exists())
                self.assertFalse(MAPPING.provenance_path(self.output).exists())
                self.assertFalse(self.raw.exists())
                self.assertEqual(Path(str(self.raw) + ".part").read_bytes(), changed + b"[success]\n")

    def test_changed_table_cannot_redeclare_its_own_snapshot_hash_in_provenance(self):
        self.acquire()
        changed = BODY.replace(b"P00001", b"P99999")
        self.output.write_bytes(changed)
        receipt = MAPPING.provenance_path(self.output)
        record = json.loads(receipt.read_text())
        record.update(bytes=len(changed), sha256=hashlib.sha256(changed).hexdigest())
        receipt.write_text(json.dumps(record))
        with self.assertRaisesRegex(ValueError, "default snapshot size/SHA256 mismatch"):
            MAPPING.verify_pinned(self.output)

    def test_explicit_manual_mapping_is_not_forced_to_the_default_snapshot(self):
        manual = self.folder / "explicit_supplied.tsv"
        changed = BODY.replace(b"P00001", b"P99999")
        manual.write_bytes(changed)
        identity = MAPPING.mapping_identity(manual)
        self.assertEqual(identity["mapping_sha256"], hashlib.sha256(changed).hexdigest())
        self.assertIsNone(identity["provenance_sha256"])

    def test_generated_name_requires_provenance_and_full_statistics(self):
        self.output.write_bytes(BODY)
        with self.assertRaises(OSError):
            MAPPING.mapping_identity(self.output)
        self.output.unlink()
        self.acquire()
        receipt = MAPPING.provenance_path(self.output)
        record = json.loads(receipt.read_text())
        record["statistics"] = {}
        receipt.write_text(json.dumps(record))
        with self.assertRaisesRegex(ValueError, "full-export statistics"):
            MAPPING.verify_pinned(self.output)

    def test_concurrent_writer_is_rejected_without_network(self):
        descriptor = os.open(self.folder, os.O_RDONLY)
        try:
            fcntl.flock(descriptor, fcntl.LOCK_EX | fcntl.LOCK_NB)
            with mock.patch.object(MAPPING, "source_metadata", side_effect=AssertionError("network")):
                with self.assertRaisesRegex(ValueError, "Another mapping writer"):
                    MAPPING.acquire(self.output, self.raw)
        finally:
            os.close(descriptor)

    def test_orphaned_receipt_published_raw_and_broken_links_are_preserved(self):
        for kind in ("receipt", "raw", "output-link", "receipt-link", "raw-link"):
            with self.subTest(kind=kind):
                target = (MAPPING.provenance_path(self.output) if kind.startswith("receipt")
                          else self.raw if kind.startswith("raw") else self.output)
                if kind.endswith("link"):
                    target.symlink_to(self.folder / "missing")
                else:
                    target.write_text("retained provenance or publisher response\n")
                with mock.patch.object(MAPPING, "source_metadata", side_effect=AssertionError("network")):
                    with self.assertRaises((ValueError, OSError)):
                        MAPPING.acquire(self.output, self.raw)
                self.assertTrue(target.exists() or target.is_symlink())
                target.unlink()

    def write_receipt(self, artifact, mapping, anchors=()):
        captured = MAPPING.mapping_identity(mapping, anchors)
        args = ["biomart_nm_to_enst.py", "receipt-write", str(artifact), str(mapping),
                "--expected-identity", json.dumps(captured)]
        for anchor in anchors:
            args.extend(("--input", str(anchor)))
        with mock.patch("sys.argv", args), contextlib.redirect_stderr(io.StringIO()):
            self.assertEqual(MAPPING.main(), 0)

    def test_receipt_detects_mapping_anchor_and_partial_output_changes(self):
        manual = self.folder / "explicit_supplied.tsv"
        manual.write_bytes(BODY)
        anchor = self.folder / "HC_spots.hg19.tsv"
        anchor.write_text("original anchors\n")
        artifact = self.folder / "motifs.pkl"
        artifact.write_bytes(b"completed scientific output")
        self.write_receipt(artifact, manual, [anchor])
        MAPPING.check_receipt(artifact, manual, [anchor])
        for path, changed in ((manual, BODY + b"ENSG000009\tENST000009\t\t\tNM_000009\n"),
                              (anchor, b"new anchors\n"), (artifact, b"partial")):
            original = path.read_bytes()
            path.write_bytes(changed)
            with self.assertRaisesRegex(ValueError, "Cannot reuse or overwrite"):
                MAPPING.check_receipt(artifact, manual, [anchor])
            path.write_bytes(original)

    def test_receipt_never_certifies_mapping_changed_during_scientific_command(self):
        manual = self.folder / "explicit_supplied.tsv"
        manual.write_bytes(BODY)
        captured = MAPPING.mapping_identity(manual)
        artifact = self.folder / "pm1.pkl"
        artifact.write_text("scientific output\n")
        manual.write_bytes(BODY.replace(b"ENST000001", b"ENST000009"))
        args = ["biomart_nm_to_enst.py", "receipt-write", str(artifact), str(manual),
                "--expected-identity", json.dumps(captured)]
        with mock.patch("sys.argv", args), contextlib.redirect_stderr(io.StringIO()):
            self.assertEqual(MAPPING.main(), 1)
        self.assertFalse(Path(str(artifact) + ".nm_to_enst.json").exists())


class GRCh38BioMartMappingTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory(prefix="biomart_grch38_", dir=os.environ.get("PRIVA_TEST_TMPDIR", os.environ.get("TMPDIR")))
        self.addCleanup(self.temp.cleanup)
        self.folder = Path(self.temp.name)
        self.output = self.folder / MAPPING.GRCH38_BASENAME
        self.raw = self.folder / "source38.raw.tsv"
        self.metadata = MAPPING.parse_metadata(REGISTRY38, DATASETS38, MAPPING.GRCH38_SOURCE_ID)
        snapshot = mock.patch.multiple(MAPPING, GRCH38_EXPECTED_SNAPSHOT_BYTES=len(BODY),
                                       GRCH38_EXPECTED_SNAPSHOT_SHA256=hashlib.sha256(BODY).hexdigest())
        snapshot.start()
        self.addCleanup(snapshot.stop)

    def acquire(self, body=BODY + b"[success]\n"):
        with mock.patch.object(MAPPING, "source_metadata", return_value=self.metadata), \
                mock.patch.object(MAPPING, "_request", return_value=io.BytesIO(body)), \
                contextlib.redirect_stdout(io.StringIO()):
            return MAPPING.acquire(self.output, self.raw, attempts=1, source_id=MAPPING.GRCH38_SOURCE_ID)

    def test_assembly_selector_has_two_fixed_sources_and_legacy_defaults(self):
        for assembly in ("hg19", "GRCh37", "GRCh37.p13"):
            selected = MAPPING.source_for_assembly(assembly)
            self.assertEqual(selected["source_id"], MAPPING.SOURCE_ID)
            self.assertEqual(selected["basename"], MAPPING.BASENAME)
            self.assertEqual(selected["relative_path"], MAPPING.RELATIVE_PATH)
            self.assertEqual((selected["mart_release"], selected["annotation_base_release"]), (116, 75))
        for assembly in ("hg38", "GRCh38", "GRCh38.p14"):
            selected = MAPPING.source_for_assembly(assembly)
            self.assertEqual(selected["source_id"], MAPPING.GRCH38_SOURCE_ID)
            self.assertEqual(selected["basename"], MAPPING.GRCH38_BASENAME)
            self.assertEqual(selected["relative_path"], MAPPING.GRCH38_RELATIVE_PATH)
            self.assertEqual((selected["mart_release"], selected["annotation_base_release"]), (113, 113))
        self.assertEqual(MAPPING.BASENAME, "ensembl_GRCh37_mart116_nm_to_enst.tsv")
        self.assertEqual(MAPPING.EXPECTED_SNAPSHOT_SHA256,
                         "f730640e2411807e7bc38dad141acafc16c3ffca5e4d8ef500976ebd586ddbc7")
        for unknown in (None, "GRCh39", "hg18"):
            with self.assertRaisesRegex(ValueError, "No explicitly pinned"):
                MAPPING.source_for_assembly(unknown)

    def test_selector_descriptors_cannot_mutate_global_source_selection(self):
        selected = MAPPING.source_for_assembly("hg38")
        selected["source_id"] = MAPPING.SOURCE_ID
        self.assertEqual(MAPPING.source_for_assembly("hg38")["source_id"], MAPPING.GRCH38_SOURCE_ID)

    def test_complete_grch38_export_preserves_order_multiplicity_and_source(self):
        record = self.acquire()
        self.assertEqual(self.output.read_bytes(), BODY)
        self.assertEqual(self.raw.read_bytes(), BODY + b"[success]\n")
        self.assertEqual(record["source_id"], MAPPING.GRCH38_SOURCE_ID)
        self.assertEqual(record["mart_release"], 113)
        self.assertEqual(record["annotation_base_release"], 113)
        self.assertEqual(record["statistics"]["duplicate_full_rows"], 1)
        self.assertEqual(record["statistics"]["nm_with_multiple_enst"], 1)
        for assembly in (None, "hg38", "GRCh38", "GRCh38.p14"):
            self.assertEqual(MAPPING.verify_pinned(self.output, assembly)["sha256"], hashlib.sha256(BODY).hexdigest())
        for assembly in ("hg19", "GRCh37", "GRCh37.p13"):
            with self.assertRaisesRegex(ValueError, "cannot be selected implicitly"):
                MAPPING.verify_pinned(self.output, assembly)

    def test_grch38_reuse_is_offline_but_never_silently_uses_default_source(self):
        self.acquire()
        with mock.patch.object(MAPPING, "source_metadata", side_effect=AssertionError("network")), \
                contextlib.redirect_stdout(io.StringIO()):
            MAPPING.acquire(self.output, self.raw, source_id=MAPPING.GRCH38_SOURCE_ID)
            with self.assertRaisesRegex(ValueError, "conflicts"):
                MAPPING.acquire(self.output, self.raw)

    def test_selected_source_cannot_publish_under_other_pinned_basename(self):
        with mock.patch.object(MAPPING, "source_metadata", side_effect=AssertionError("network")):
            with self.assertRaisesRegex(ValueError, "conflicts"):
                MAPPING.acquire(self.folder / MAPPING.BASENAME, self.raw, source_id=MAPPING.GRCH38_SOURCE_ID)
        self.assertFalse(self.raw.exists())

    def test_known_filename_cannot_disagree_with_provenance_source(self):
        self.acquire()
        other = self.folder / MAPPING.BASENAME
        self.output.rename(other)
        MAPPING.provenance_path(self.output).rename(MAPPING.provenance_path(other))
        with self.assertRaisesRegex(ValueError, "filename conflicts"):
            MAPPING.verify_pinned(other)

    def test_grch38_release_rollover_and_wrong_assembly_metadata_are_rejected(self):
        with self.assertRaisesRegex(ValueError, "source changed"):
            MAPPING.parse_metadata(REGISTRY38.replace("113", "114"), DATASETS38, MAPPING.GRCH38_SOURCE_ID)
        with self.assertRaisesRegex(ValueError, "dataset/assembly changed"):
            MAPPING.parse_metadata(REGISTRY38, DATASETS, MAPPING.GRCH38_SOURCE_ID)
        with self.assertRaisesRegex(ValueError, "source changed"):
            MAPPING.parse_metadata(REGISTRY38, DATASETS38)

    def test_grch38_changed_order_or_values_cannot_publish_a_new_snapshot(self):
        lines = BODY.splitlines(keepends=True)
        reordered = b"".join([lines[0], lines[3], lines[1], lines[2], *lines[4:]])
        for changed in (reordered, BODY.replace(b"P00001", b"P99999")):
            with self.assertRaisesRegex(ValueError, "default snapshot size/SHA256 mismatch"):
                self.acquire(changed + b"[success]\n")
            self.assertFalse(self.output.exists())
            self.assertFalse(self.raw.exists())

    def test_grch38_table_cannot_redeclare_its_own_hash_in_provenance(self):
        self.acquire()
        self.output.write_bytes(BODY.replace(b"P00001", b"P99999"))
        path = MAPPING.provenance_path(self.output)
        record = json.loads(path.read_text())
        record.update(bytes=self.output.stat().st_size, sha256=MAPPING.sha256(self.output))
        path.write_text(json.dumps(record))
        with self.assertRaisesRegex(ValueError, "default snapshot size/SHA256 mismatch"):
            MAPPING.verify_pinned(self.output)

    def test_unknown_provenance_source_is_never_a_fallback(self):
        self.acquire()
        path = MAPPING.provenance_path(self.output)
        record = json.loads(path.read_text()); record["source_id"] = "ensembl-grch38-current"
        path.write_text(json.dumps(record))
        with self.assertRaisesRegex(ValueError, "Unknown explicitly pinned"):
            MAPPING.verify_pinned(self.output)

    def test_generated_grch38_name_requires_provenance_for_derived_identity(self):
        self.output.write_bytes(BODY)
        with self.assertRaises(OSError):
            MAPPING.mapping_identity(self.output)

    def test_fetch_cli_dispatches_the_explicit_grch38_source(self):
        args = ["biomart_nm_to_enst.py", "fetch", "--source", MAPPING.GRCH38_SOURCE_ID,
                "--output", str(self.output), "--raw-output", str(self.raw)]
        with mock.patch("sys.argv", args), mock.patch.object(MAPPING, "acquire") as acquire:
            self.assertEqual(MAPPING.main(), 0)
        self.assertEqual(acquire.call_args.kwargs["source_id"], MAPPING.GRCH38_SOURCE_ID)

    def test_request_uses_selected_archive_and_rejects_cross_host_redirect(self):
        response = mock.Mock(url="https://grch37.ensembl.org/biomart/martservice")
        with mock.patch.object(MAPPING, "urlopen", return_value=response) as request:
            with self.assertRaisesRegex(ValueError, "redirected away"):
                MAPPING._request({"type": "registry"}, 1, MAPPING.GRCH38_SOURCE_ID)
        self.assertTrue(request.call_args.args[0].full_url.startswith(MAPPING.GRCH38_ENDPOINT + "?"))
        response.close.assert_called_once()


if __name__ == "__main__":
    unittest.main()
