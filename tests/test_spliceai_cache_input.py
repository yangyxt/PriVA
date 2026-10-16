"""Synthetic header/routing checks; no score conversion or full-cache claims."""
import gzip
import hashlib
import importlib.util
import json
import os
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest
from unittest import mock

ROOT = Path(__file__).resolve().parents[1]
SPEC = importlib.util.spec_from_file_location("spliceai_cache_input", ROOT / "scripts/spliceai_cache_input.py")
CACHE = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(CACHE)
VCF = "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
NATIVE = "\t".join(CACHE.native_columns()) + "\n"


class Fixture(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory(prefix="spliceai_cache_input_")
        self.addCleanup(self.temp.cleanup)
        self.folder = Path(self.temp.name)
        self.snv = self.write("selected-snv", VCF)
        self.indel = self.write("selected-indel", VCF)
        self.config = self.write("config.yaml", "# untouched\nassembly: hg38\nother: [one, two]\n")

    def write(self, name, text, compressed=False):
        path = self.folder / name
        value = text.encode() if isinstance(text, str) else text
        path.write_bytes(gzip.compress(value) if compressed else value)
        return path

    def header(self, metadata="", columns=None):
        return "##fileformat=VCFv4.2\n" + metadata + (columns or "\t".join(CACHE.VCF_COLUMNS)) + "\n"

    def full_manifest(self):
        files = []
        for kind, chrom in sorted(CACHE.PUBLISHED):
            path = self.write(f"SpliceAI2_GRCh38_{kind}_{chrom}.tsv.gz", NATIVE, compressed=True)
            files.append({"kind": kind, "chrom": chrom, "path": str(path), "bytes": path.stat().st_size,
                          "sha256": hashlib.sha256(path.read_bytes()).hexdigest()})
        self.snv = self.folder / "SpliceAI2_GRCh38_snv_chr1.tsv.gz"
        self.indel = self.folder / "SpliceAI2_GRCh38_indel_chr1.tsv.gz"
        self.manifest = self.folder / CACHE.MANIFEST_NAME
        self.document = {"schema_version": 1, "resource": "SpliceAI2", "version": "2.0", "assembly": "GRCh38",
                         "status": "downloads_complete", "confirmed_by": "synthetic-fixture",
                         "source_revision": CACHE.REVISION, "files": files}
        self.save_manifest()

    def save_manifest(self):
        self.manifest.write_text(json.dumps(self.document))


class HeaderTests(Fixture):
    def test_plain_vcf_with_gzip_suffix_and_gzip_without_suffix(self):
        for name, compressed in (("plain.vcf.gz", False), ("compressed.anything", True)):
            with self.subTest(name=name):
                path = self.write(name, VCF, compressed)
                self.assertEqual(CACHE.inspect_header(path), {"format": "vcf_v1", "assembly": None})

    def test_native_exact_header_hash_and_gzip_magic(self):
        self.assertEqual(len(CACHE.native_columns()), 84)
        self.assertEqual(hashlib.sha256(NATIVE.encode()).hexdigest(),
                         "1feb799bed12dae1030e64b03930eab470340848aa4a7e0166072a140e7b0d66")
        for compressed in (False, True):
            path = self.write("table", NATIVE, compressed)
            self.assertEqual(CACHE.inspect_header(path), {"format": "spliceai2_tsv", "assembly": "GRCh38"})

    def test_native_reordered_missing_extra_or_comment_header_rejected(self):
        columns = CACHE.native_columns()
        variants = [columns[:-1], columns + ["extra"], [columns[1], columns[0]] + columns[2:]]
        for fields in variants:
            with self.subTest(fields=fields[:2], length=len(fields)):
                with self.assertRaisesRegex(ValueError, "Unsupported cache header"):
                    CACHE.inspect_header(self.write("bad", "\t".join(fields) + "\n"))
        with self.assertRaises(ValueError):
            CACHE.inspect_header(self.write("bad", "#" + NATIVE))

    def test_supported_vcf_samples_and_crlf(self):
        text = self.header('##INFO=<ID=SpliceAI,Number=.,Type=String,Description="a, b">\n',
                           "\t".join(CACHE.VCF_COLUMNS + ["FORMAT", "SAMPLE_1", "SAMPLE_2"]))
        self.assertEqual(CACHE.inspect_header(self.write("samples", text.replace("\n", "\r\n")))["format"], "vcf_v1")

    def test_bad_vcf_columns_or_header_order_rejected(self):
        cases = [VCF.replace("#CHROM\tPOS", "#CHROM\tPOSITION"),
                 self.header(columns="\t".join(CACHE.VCF_COLUMNS + ["FORMAT"])),
                 self.header(columns="\t".join(CACHE.VCF_COLUMNS + ["FORMAT", "same", "same"])),
                 VCF.split("\n", 1)[1], self.header("##fileformat=VCFv4.2\n"),
                 "##fileformat=VCFv4.2\n1\t10\t.\tA\tC\t.\t.\t.\n"]
        for text in cases:
            with self.subTest(text=text[:90]), self.assertRaises(ValueError):
                CACHE.inspect_header(self.write("bad", text))

    def test_truncated_or_invalid_metadata_rejected(self):
        for metadata in ('##INFO=<ID=x\n', '##INFO=<ID=x,Description="unclosed>\n', '##broken\n'):
            with self.subTest(metadata=metadata), self.assertRaises(ValueError):
                CACHE.inspect_header(self.write("bad", self.header(metadata)))

    def test_assembly_reference_and_contig_conflicts(self):
        for metadata, assembly in (("##reference=file:///refs/hg19.fa\n", "GRCh37"),
                                   ('##contig=<ID=1,length=10,assembly="GRCh38.p14">\n', "GRCh38")):
            self.assertEqual(CACHE.inspect_header(self.write("metadata", self.header(metadata)))["assembly"], assembly)
        with self.assertRaisesRegex(ValueError, "Conflicting assembly"):
            CACHE.inspect_header(self.write("conflict", self.header("##reference=hg19\n##assembly=GRCh38\n")))

    def test_projection_markers_and_version_conflicts(self):
        for metadata in ("##source=PriVA_SpliceAI2_projection\n", "##SpliceAI2Version=2\n",
                         "##SpliceAIVersion=2.0\n"):
            self.assertEqual(CACHE.inspect_header(self.write("v2", self.header(metadata)))["format"], "vcf_spliceai2")
        for metadata in ("##SpliceAI2Version=1\n", "##SpliceAI2Version=2\n##SpliceAIVersion=1\n",
                         "##SpliceAI2Version=2\n##assembly=hg19\n"):
            with self.subTest(metadata=metadata), self.assertRaises(ValueError):
                CACHE.inspect_header(self.write("conflict", self.header(metadata)))

    def test_corrupt_gzip_and_missing_header_newlines(self):
        cases = [b"\x1f\x8bnot-a-valid-gzip", gzip.compress(VCF.encode())[:-6],
                 VCF.rstrip("\n").encode(), gzip.compress(NATIVE.rstrip("\n").encode()),
                 b"##fileformat=VCFv4.2\n##source=\xff\n", VCF.replace("POS", "PO\x00S").encode()]
        for value in cases:
            with self.subTest(prefix=value[:12]), self.assertRaises((ValueError, OSError)):
                CACHE.inspect_header(self.write("bad", value))

    def test_partial_input_and_partial_symlink_target_rejected(self):
        partial = self.write("source.partial.gz", VCF, True)
        alias = self.folder / "alias"
        alias.symlink_to(partial)
        for path in (partial, alias):
            with self.subTest(path=path.name), self.assertRaisesRegex(ValueError, "Incomplete"):
                CACHE.inspect_header(path)

    def test_bounded_line_bytes_total_bytes_and_line_count(self):
        for name, limit, text in (("MAX_LINE", 10, VCF),
                                  ("MAX_HEADER", 50, VCF),
                                  ("MAX_LINES", 2, self.header("##source=test\n"))):
            with self.subTest(limit=name), mock.patch.object(CACHE, name, limit), self.assertRaisesRegex(ValueError, "bounded"):
                CACHE.inspect_header(self.write("bounded", text))

    def test_later_body_is_not_scanned_or_claimed_validated(self):
        # This malformed body is deliberately outside the header inspector's
        # scope: neither score interpretation nor complete gzip CRC is claimed.
        path = self.write("body", VCF + "NOT_A_VARIANT\n" + "x" * 100000)
        self.assertEqual(CACHE.inspect_header(path)["format"], "vcf_v1")


class PairAndInstallTests(Fixture):
    def test_existing_v1_vcf_both_assemblies_pass_without_conversion(self):
        for assembly in ("hg19", "GRCh37", "hg38", "GRCh38"):
            self.assertEqual(CACHE.validate_pair(assembly, self.snv, self.indel)["format"], "vcf_v1")

    def test_mixed_pair_and_explicit_assembly_conflict_fail(self):
        self.indel.write_text(NATIVE)
        with self.assertRaisesRegex(ValueError, "Mixed"):
            CACHE.validate_pair("hg38", self.snv, self.indel)
        self.indel.write_text(self.header("##SpliceAI2Version=2\n"))
        with self.assertRaisesRegex(ValueError, "Mixed"):
            CACHE.validate_pair("hg38", self.snv, self.indel)
        self.snv.write_text(self.indel.read_text())
        with self.assertRaisesRegex(ValueError, "assembly conflicts"):
            CACHE.validate_pair("hg19", self.snv, self.indel)
        self.snv.write_text(self.header("##reference=hg19\n"))
        self.indel.write_text(VCF)
        with self.assertRaisesRegex(ValueError, "assembly conflicts"):
            CACHE.validate_pair("hg38", self.snv, self.indel)

    def test_vcf_registration_preserves_unrelated_text_comments_and_mode(self):
        original = '# before\r\nassembly: hg38\r\nspliceai_snv_prescore: "old" # snv comment\r\nother: [one, two] # unchanged\r\nspliceai_indel_prescore: old2 # indel comment\r\n'
        self.config.write_bytes(original.encode())
        self.config.chmod(0o640)
        with mock.patch.object(CACHE.subprocess, "run") as run:
            CACHE.install(self.config, self.folder, self.snv, self.indel)
            run.assert_not_called()
        expected = original.replace('"old"', json.dumps(str(self.snv))).replace("old2", json.dumps(str(self.indel)))
        self.assertEqual(self.config.read_bytes(), expected.encode())
        self.assertEqual(self.config.stat().st_mode & 0o777, 0o640)

    def test_existing_same_values_are_byte_and_identity_noop(self):
        text = f"assembly: hg38\nspliceai_snv_prescore: '{self.snv}' # keep\nspliceai_indel_prescore: {self.indel}\n"
        self.config.write_text(text)
        before = CACHE._identity(self.config)
        CACHE.install(self.config, self.folder, self.snv, self.indel)
        self.assertEqual(self.config.read_text(), text)
        self.assertEqual(CACHE._identity(self.config), before)

    def test_missing_keys_added_and_empty_scalar_keys_filled(self):
        for existing in ("", "spliceai_snv_prescore: # keep\nspliceai_indel_prescore:\n"):
            self.config.write_text("assembly: hg38\n" + existing)
            CACHE.install(self.config, self.folder, self.snv, self.indel)
            import yaml
            data = yaml.safe_load(self.config.read_text())
            self.assertEqual(data["spliceai_snv_prescore"], str(self.snv))
            self.assertEqual(data["spliceai_indel_prescore"], str(self.indel))
            if existing:
                self.assertIn("# keep", self.config.read_text())

    def test_bad_second_header_keeps_config_and_never_dispatches(self):
        self.indel.write_text("bad\n")
        before = self.config.read_bytes()
        with mock.patch.object(CACHE.subprocess, "run") as run, self.assertRaises(ValueError):
            CACHE.install(self.config, self.folder, self.snv, self.indel)
        self.assertEqual(self.config.read_bytes(), before)
        run.assert_not_called()

    def test_config_assembly_change_during_validation_is_rejected(self):
        original = CACHE.validate_pair

        def mutate(*args):
            result = original(*args)
            self.config.write_text("assembly: hg19\n")
            return result

        with mock.patch.object(CACHE, "validate_pair", side_effect=mutate), self.assertRaisesRegex(ValueError, "Config changed"):
            CACHE.install(self.config, self.folder, self.snv, self.indel)
        self.assertEqual(self.config.read_text(), "assembly: hg19\n")

    def test_cache_change_during_pair_validation_is_rejected(self):
        original = CACHE.validate_pair

        def mutate(*args):
            result = original(*args)
            self.snv.write_text(self.header("##reference=hg19\n"))
            return result

        before = self.config.read_bytes()
        with mock.patch.object(CACHE, "validate_pair", side_effect=mutate), self.assertRaisesRegex(ValueError, "Cache input changed"):
            CACHE.install(self.config, self.folder, self.snv, self.indel)
        self.assertEqual(self.config.read_bytes(), before)

    def test_cli_malformed_yaml_is_clean_failure(self):
        self.config.write_text("assembly: [unterminated\n")
        result = subprocess.run([sys.executable, str(ROOT / "scripts/spliceai_cache_input.py"), "install",
                                 "--config", str(self.config), "--cache-dir", str(self.folder),
                                 "--snv", str(self.snv), "--indel", str(self.indel)], capture_output=True, text=True)
        self.assertEqual(result.returncode, 2)
        self.assertIn("Invalid YAML config", result.stderr)
        self.assertNotIn("Traceback", result.stderr)


class ManifestTests(Fixture):
    def setUp(self):
        super().setUp()
        self.full_manifest()

    def test_full47_membership_and_symlink_selected_paths(self):
        self.assertEqual(CACHE.discover_manifest(self.snv, self.indel), self.manifest)
        alias = self.folder / "snv-alias"
        alias.symlink_to(self.snv)
        self.assertEqual(CACHE.discover_manifest(alias, self.indel), self.manifest)

    def test_native_dispatch_exact_arguments_no_direct_config_mutation(self):
        before = self.config.read_bytes()
        with mock.patch.object(CACHE.subprocess, "run") as run:
            result = CACHE.install(self.config, self.folder / "cache", self.snv, self.indel)
        self.assertEqual(result["format"], "spliceai2_tsv")
        run.assert_called_once_with(["bash", str(ROOT / "scripts/install_utils.sh"), "SpliceAI2_install", str(self.config),
                                     str(self.manifest), str(self.folder / "cache/SpliceAI2/GRCh38"), "priva"], check=True)
        self.assertEqual(self.config.read_bytes(), before)

    def test_absent_manifest_no_conversion_or_fallback(self):
        self.manifest.unlink()
        with mock.patch.object(CACHE.subprocess, "run") as run, self.assertRaises(OSError):
            CACHE.install(self.config, self.folder, self.snv, self.indel)
        run.assert_not_called()

    def test_wrong_role_or_unlisted_native_input_rejected(self):
        with self.assertRaisesRegex(ValueError, "manifest role"):
            CACHE.discover_manifest(self.indel, self.snv)
        extra = self.write("unlisted", NATIVE)
        with self.assertRaisesRegex(ValueError, "manifest role"):
            CACHE.discover_manifest(extra, self.indel)

    def test_subset_duplicate_role_and_unpublished_indel_y_rejected(self):
        original = json.loads(json.dumps(self.document))
        variants = [original["files"][:-1], original["files"][:-1] + [original["files"][0]],
                    [{**item, "chrom": "chrY"} if item["kind"] == "indel" and item["chrom"] == "chr1" else item
                     for item in original["files"]]]
        for files in variants:
            self.document["files"] = files
            self.save_manifest()
            with self.subTest(count=len(files)), self.assertRaises(ValueError):
                CACHE.discover_manifest(self.snv, self.indel)

    def test_unconfirmed_or_wrong_source_identity_rejected(self):
        original = dict(self.document)
        for key, value in (("confirmed_by", " "), ("status", "downloading"), ("assembly", "GRCh37"),
                           ("source_revision", "different"), ("schema_version", True), ("version", "1.3")):
            self.document = {**original, key: value}
            self.save_manifest()
            with self.subTest(key=key), self.assertRaises(ValueError):
                CACHE.discover_manifest(self.snv, self.indel)

    def test_file_size_sha_shape_and_incomplete_paths_rejected(self):
        original = dict(self.document["files"][0])
        for key, value in (("bytes", original["bytes"] + 1), ("bytes", True), ("sha256", "missing"),
                           ("path", "relative.tsv.gz"), ("kind", [])):
            self.document["files"][0] = {**original, key: value}
            self.save_manifest()
            with self.subTest(key=key, value=value), self.assertRaises(ValueError):
                CACHE.discover_manifest(self.snv, self.indel)

    def test_duplicate_json_key_rejected(self):
        text = self.manifest.read_text()
        self.manifest.write_text(text[:-1] + ', "status":"downloads_complete"}')
        with self.assertRaisesRegex(ValueError, "Duplicate manifest key"):
            CACHE.discover_manifest(self.snv, self.indel)

    def test_native_hg19_never_dispatches(self):
        self.config.write_text("assembly: hg19\n")
        with mock.patch.object(CACHE.subprocess, "run") as run, self.assertRaisesRegex(ValueError, "assembly conflicts"):
            CACHE.install(self.config, self.folder, self.snv, self.indel)
        run.assert_not_called()

    def test_native_config_mutation_before_dispatch_rejected(self):
        original = CACHE.discover_manifest

        def mutate(*args):
            result = original(*args)
            self.config.write_text("assembly: hg19\n")
            return result

        with mock.patch.object(CACHE, "discover_manifest", side_effect=mutate), mock.patch.object(CACHE.subprocess, "run") as run:
            with self.assertRaisesRegex(ValueError, "Config changed"):
                CACHE.install(self.config, self.folder, self.snv, self.indel)
        run.assert_not_called()


if __name__ == "__main__":
    unittest.main()
