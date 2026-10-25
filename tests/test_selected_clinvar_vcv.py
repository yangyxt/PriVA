"""Explicit VCV acquisition/identity only; no scientific filtering or network tests."""
from contextlib import redirect_stdout
import gzip
import hashlib
import io
import json
import os
from pathlib import Path
import sys
import tempfile
import unittest
from unittest.mock import patch


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
import build_gene_nonlof_mechanism_cache as builder
import inject_clinvar_into_gofcards as injector


VCV = (b"<ClinVarVariationRelease ReleaseDate='2026-10-11'>"
       b"<VariationArchive VariationID='1' RecordType='classified'><RecordStatus>current</RecordStatus></VariationArchive>"
       b"<VariationArchive VariationID='2' RecordType='included'><RecordStatus>replaced</RecordStatus></VariationArchive>"
       b"</ClinVarVariationRelease>")


class SelectedClinVarVCV(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory(
            prefix="priva_selected_vcv_", dir=os.environ.get("PRIVA_TEST_TMPDIR", os.environ.get("TMPDIR")))
        self.addCleanup(self.temporary.cleanup)
        self.folder = Path(self.temporary.name)
        self.selected = self.folder / "selected source's.xml.gz"
        with gzip.open(self.selected, "wb") as stream:
            stream.write(VCV)

    def args(self, selected=True, fetch_sources_only=True):
        arguments = ["builder", "--cache-dir", str(self.folder / "cache"),
                     "--shared-raw-dir", str(self.folder / "raw")]
        if selected:
            arguments += ["--clinvar-vcv", str(self.selected)]
        if fetch_sources_only:
            arguments.append("--fetch-sources-only")
        with patch.object(sys, "argv", arguments):
            return builder.parse_args()

    def test_genuine_gzip_and_plain_xml_identity_without_publisher_claim(self):
        plain = self.folder / "selected.xml"
        plain.write_bytes(VCV)
        for path in (self.selected, plain):
            with self.subTest(path=path):
                before = path.read_bytes()
                metadata = builder.selected_clinvar_vcv_metadata(path)
                self.assertEqual(path.read_bytes(), before)
                self.assertEqual(metadata["raw_path"], str(path))
                self.assertEqual(metadata["sha256"], hashlib.sha256(before).hexdigest())
                self.assertEqual(metadata["md5"], hashlib.md5(before).hexdigest())
                self.assertEqual(metadata["size_bytes"], len(before))
                self.assertEqual(metadata["variation_archive_records"], 2)
                self.assertEqual(metadata["declared_release_date"], "2026-10-11")
                self.assertIsNone(metadata["url"])
                self.assertIsNone(metadata["remote_md5"])
                self.assertFalse(metadata["publisher_identity_verified"])
                self.assertEqual(metadata["source_mode"], "explicit_local_input")

    def test_invalid_xml_root_empty_records_and_missing_variation_id_are_retained(self):
        for content in (b"not XML", b"<interprodb/>", b"<ClinVarVariationRelease/>",
                        b"<ClinVarVariationRelease><VariationArchive/></ClinVarVariationRelease>",
                        VCV[:-5], VCV + b"trailing corruption"):
            with self.subTest(content=content[:50]):
                with gzip.open(self.selected, "wb") as stream:
                    stream.write(content)
                before = self.selected.read_bytes()
                with self.assertRaises((ValueError, builder.ET.ParseError)):
                    builder.selected_clinvar_vcv_metadata(self.selected)
                self.assertEqual(self.selected.read_bytes(), before)

    def test_truncated_gzip_fails_instead_of_accepting_an_xml_prefix(self):
        self.selected.write_bytes(self.selected.read_bytes()[:-5])
        before = self.selected.read_bytes()
        with self.assertRaises((EOFError, OSError)):
            builder.selected_clinvar_vcv_metadata(self.selected)
        self.assertEqual(self.selected.read_bytes(), before)

    def test_changed_during_validation_fails_identity_check(self):
        real_parse = builder.ET.iterparse

        def changed(*arguments, **kwargs):
            yield from real_parse(*arguments, **kwargs)
            stat = self.selected.stat()
            os.utime(self.selected, ns=(stat.st_atime_ns, stat.st_mtime_ns + 1_000_000))

        with patch.object(builder.ET, "iterparse", changed):
            with self.assertRaisesRegex(ValueError, "changed during validation"):
                builder.selected_clinvar_vcv_metadata(self.selected)

    def test_missing_or_corrupt_selected_input_fails_before_all_acquisition_and_cache_writes(self):
        for missing in (True, False):
            with self.subTest(missing=missing):
                if missing:
                    self.selected.unlink(missing_ok=True)
                else:
                    self.selected.write_bytes(b"retained corrupt input")
                with patch.object(builder, "parse_args", return_value=self.args()), \
                     patch.object(builder, "fetch_static_source") as static, \
                     patch.object(builder, "fetch_panelapp") as panel, \
                     patch.object(builder, "fetch_clinvar_vcv") as weekly:
                    with self.assertRaises((ValueError, OSError)):
                        builder.main()
                    static.assert_not_called()
                    panel.assert_not_called()
                    weekly.assert_not_called()
                self.assertFalse((self.folder / "cache").exists())
                self.assertFalse((self.folder / "raw").exists())

    def test_sources_only_records_exact_custom_identity_and_skips_weekly_fetch(self):
        metadata = {"sha256": "other-source-fixture", "status": "stubbed"}
        with patch.object(builder, "parse_args", return_value=self.args()), \
             patch.object(builder, "fetch_static_source", return_value=metadata), \
             patch.object(builder, "fetch_panelapp", return_value=metadata), \
             patch.object(builder, "fetch_clinvar_vcv") as weekly, \
             patch.object(builder, "parse_all_sources") as scientific_parser, redirect_stdout(io.StringIO()):
            self.assertEqual(builder.main(), 0)
            weekly.assert_not_called()
            scientific_parser.assert_not_called()
        manifest = json.loads((self.folder / "cache" / builder.NONLOF_SOURCE_MANIFEST_FILENAME).read_text())
        selected = manifest["sources"][builder.CLINVAR_VCV_SPEC.name]
        self.assertEqual(selected["raw_path"], str(self.selected))
        self.assertEqual(selected["sha256"], hashlib.sha256(self.selected.read_bytes()).hexdigest())
        self.assertEqual(selected["variation_archive_records"], 2)
        self.assertEqual(manifest["selected_clinvar_vcv_path"], str(self.selected))
        self.assertNotIn(builder.CLINVAR_VCV_SPEC.raw_filename, manifest["expected_raw_filenames"])
        self.assertFalse((self.folder / "raw" / builder.CLINVAR_VCV_SPEC.raw_filename).exists())

    def test_omitted_optional_path_preserves_original_weekly_fetch_contract(self):
        metadata = {"sha256": "stub-source-identity", "status": "stubbed"}
        with patch.object(builder, "parse_args", return_value=self.args(selected=False)), \
             patch.object(builder, "fetch_static_source", return_value=metadata), \
             patch.object(builder, "fetch_panelapp", return_value=metadata), \
             patch.object(builder, "fetch_clinvar_vcv", return_value=metadata) as weekly, \
             patch.object(builder, "selected_clinvar_vcv_metadata") as selected, redirect_stdout(io.StringIO()):
            self.assertEqual(builder.main(), 0)
            selected.assert_not_called()
            weekly.assert_called_once()
            self.assertEqual(weekly.call_args.kwargs["raw_dir"], self.folder / "raw")

    def test_full_builder_records_selected_path_without_changing_parser_calls(self):
        args = self.args(fetch_sources_only=False)
        args.gofcards_exact_variants = self.folder / "gofcards.json.gz"
        args.hgnc_table = self.folder / "hgnc.tsv"
        args.output_schema = self.folder / "schema.json"
        for path in (args.gofcards_exact_variants, args.hgnc_table, args.output_schema):
            path.write_text("operational identity fixture\n")
        metadata = {"sha256": "stub-source-identity", "status": "stubbed"}
        canonical = {"_meta": {"total_genes": 0, "unmapped_symbols": []}}
        with patch.object(builder, "parse_args", return_value=args), \
             patch.object(builder, "load_gofcards_cache", return_value={"metadata": {"clinvar": {
                 "source": str(self.selected), "source_sha256": hashlib.sha256(self.selected.read_bytes()).hexdigest()}}}), \
             patch.object(builder, "fetch_static_source", return_value=metadata), \
             patch.object(builder, "fetch_panelapp", return_value=metadata), \
             patch.object(builder, "fetch_clinvar_vcv") as weekly, \
             patch.object(builder, "parse_all_sources", return_value="all-original-source-rows") as parser, \
             patch.object(builder, "load_hgnc_mapping", return_value={}), \
             patch.object(builder, "load_gofcards_exact_records", return_value=({}, {})), \
             patch.object(builder, "build_nonlof_assertions_json", return_value=canonical) as build, \
             patch.object(builder, "validate_canonical_json"), redirect_stdout(io.StringIO()):
            self.assertEqual(builder.main(), 0)
            weekly.assert_not_called()
            parser.assert_called_once_with(raw_dir=self.folder / "raw")
            build.assert_called_once_with("all-original-source-rows", {}, gofcards_exact_by_variant={})
        manifest = json.loads((self.folder / "cache" / builder.NONLOF_SOURCE_MANIFEST_FILENAME).read_text())
        self.assertEqual(manifest["build_inputs"]["selected_clinvar_vcv_path"], str(self.selected))
        self.assertEqual(manifest["build_inputs"]["source_sha256"][builder.CLINVAR_VCV_SPEC.name],
                         hashlib.sha256(self.selected.read_bytes()).hexdigest())

    def test_explicit_full_build_rejects_injected_xml_path_or_digest_mismatch_before_acquisition(self):
        for source, digest in ((str(self.folder / "different.xml.gz"), hashlib.sha256(self.selected.read_bytes()).hexdigest()),
                               (str(self.selected), "0" * 64), (str(self.selected), None)):
            with self.subTest(source=source, digest=digest):
                with patch.object(builder, "parse_args", return_value=self.args(fetch_sources_only=False)), \
                     patch.object(builder, "load_gofcards_cache", return_value={"metadata": {"clinvar": {"source": source, "source_sha256": digest}}}), \
                     patch.object(builder, "fetch_static_source") as acquire, \
                     patch.object(builder, "parse_all_sources") as scientific:
                    with self.assertRaisesRegex(ValueError, "does not match the injected"):
                        builder.main()
                    acquire.assert_not_called()
                    scientific.assert_not_called()
                self.assertFalse((self.folder / "cache").exists())

    def test_same_path_changed_xml_cannot_relabel_existing_injected_evidence(self):
        original_sha = hashlib.sha256(self.selected.read_bytes()).hexdigest()
        with gzip.open(self.selected, "wb") as handle:
            handle.write(VCV.replace(b"2026-10-11", b"2026-10-12"))
        current = builder.selected_clinvar_vcv_metadata(self.selected)
        self.assertNotEqual(current["sha256"], original_sha)
        with patch.object(builder, "load_gofcards_cache", return_value={"metadata": {"clinvar": {
            "source": str(self.selected), "source_sha256": original_sha}}}):
            with self.assertRaisesRegex(ValueError, "does not match the injected"):
                builder.require_selected_clinvar_injection(current, self.folder / "injected.json.gz")


class CustomVCVInjection(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory(
            prefix="priva_vcv_injection_", dir=os.environ.get("PRIVA_TEST_TMPDIR", os.environ.get("TMPDIR")))
        self.addCleanup(self.temporary.cleanup)
        self.folder = Path(self.temporary.name)
        self.selected = self.folder / "selected.xml.gz"
        with gzip.open(self.selected, "wb") as stream:
            stream.write(VCV)
        self.normalized = self.folder / "normalized.json.gz"
        self.normalized.write_text("normalized cache fixture\n")
        self.output = self.folder / "injected.json.gz"
        self.output.write_bytes(b"previous published output retained\n")

    def call(self, digest=None, parser=None):
        arguments = ["--gofcards-cache", str(self.normalized), "--clinvar-vcv", str(self.selected),
                     "--out-json", str(self.output)]
        if digest is not None:
            arguments += ["--clinvar-vcv-sha256", digest]
        with patch.object(injector, "load_gofcards_cache", return_value={"metadata": {}}), \
             patch.object(injector, "index_gofcards_variants", return_value={}), \
             patch.object(injector, "stream_parse_clinvar_vcv", side_effect=parser or (lambda *_args, **_kwargs: ([], {}))) as stream, \
             patch.object(injector, "inject", return_value={}):
            result = injector.main(arguments)
        return result, stream

    def test_exact_custom_digest_is_recorded_and_original_parser_arguments_are_unchanged(self):
        before = self.selected.read_bytes()
        digest = hashlib.sha256(before).hexdigest()
        result, stream = self.call(digest)
        self.assertEqual(result, 0)
        stream.assert_called_once_with(self.selected, {}, min_review_stars=2, variation_id_lookup={})
        with gzip.open(self.output, "rt") as handle:
            metadata = json.load(handle)["metadata"]["clinvar"]
        self.assertEqual(metadata["source"], str(self.selected))
        self.assertEqual(metadata["source_sha256"], digest)
        self.assertEqual(self.selected.read_bytes(), before)

    def test_mismatched_custom_digest_fails_before_parser_and_keeps_published_output(self):
        before = self.output.read_bytes()
        with patch.object(injector, "stream_parse_clinvar_vcv") as stream:
            with self.assertRaisesRegex(SystemExit, "mismatch before injection"):
                self.call("0" * 64)
            stream.assert_not_called()
        self.assertEqual(self.output.read_bytes(), before)

    def test_input_mutation_during_injection_prevents_publication(self):
        before = self.output.read_bytes()
        digest = hashlib.sha256(self.selected.read_bytes()).hexdigest()

        def mutate(*_args, **_kwargs):
            self.selected.write_bytes(self.selected.read_bytes() + b"changed during parser\n")
            return [], {}

        with self.assertRaisesRegex(SystemExit, "changed during injection"):
            self.call(digest, mutate)
        self.assertEqual(self.output.read_bytes(), before)

    def test_default_weekly_mode_adds_no_custom_xml_hash_or_metadata_claim(self):
        original_sha = injector.sha256_file
        with patch.object(injector, "sha256_file", wraps=original_sha) as sha:
            result, _ = self.call()
        self.assertEqual(result, 0)
        self.assertNotIn(self.selected, [call.args[0] for call in sha.call_args_list])
        with gzip.open(self.output, "rt") as handle:
            metadata = json.load(handle)["metadata"]["clinvar"]
        self.assertNotIn("source_sha256", metadata)


if __name__ == "__main__":
    unittest.main()
