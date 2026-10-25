"""Standalone installer conversion, compiler isolation and production-reference parity."""
import importlib.util
import json
import os
from pathlib import Path
import shutil
import subprocess
import tempfile
import unittest
from unittest import mock


ROOT = Path(__file__).resolve().parents[1]
SPEC = importlib.util.spec_from_file_location("spliceai2_install_driver", ROOT / "scripts/spliceai2_to_vep_vcf.py")
DRIVER = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(DRIVER)


def native_columns():
    columns = ["chrom", "pos", "ref", "alt", "strand"]
    for effect in ("donor_gain", "donor_loss", "acceptor_gain", "acceptor_loss", "jxn_gain", "jxn_loss"):
        fields = ["delta_score", "ref_score", "alt_score"]
        fields += ["donor_dist", "acceptor_dist"] if effect.startswith("jxn") else ["dist"]
        columns.extend(f"{effect}_{field}_{rank}" for field in fields for rank in range(3))
    return columns + ["spliceai2_summary_score"]


def row(*, pos="11", ref="A", alt="C", strand="+", **changes):
    columns = native_columns()
    values = {name: "0" if "_dist_" in name else "0.00" for name in columns}
    values.update(chrom="chr1", pos=pos, ref=ref, alt=alt, strand=strand,
                  acceptor_gain_delta_score_0="0.2500", donor_gain_delta_score_0="-0.0",
                  acceptor_gain_dist_0="-3", spliceai2_summary_score="0.25")
    values.update(changes)
    return "\t".join(values[name] for name in columns)


class CompilerTests(unittest.TestCase):
    def setUp(self):
        temporary = tempfile.TemporaryDirectory(prefix="spliceai2_compiler_",
                                                dir=os.environ.get("PRIVA_TEST_TMPDIR", ROOT / "tests"))
        self.addCleanup(temporary.cleanup)
        self.folder = Path(temporary.name)

    def test_compatible_host_compiler_needs_no_package_manager(self):
        compiler = dict(path="/host/bin/rustc", version="rustc 1.90.0 (fixture)", sha256="f" * 64)
        with mock.patch.object(DRIVER, "compiler_identity", return_value=compiler), \
                mock.patch.object(DRIVER.subprocess, "run") as run:
            self.assertEqual(DRIVER.installation_compiler(self.folder), compiler)
        run.assert_not_called()

    def test_compiler_identity_distinguishes_shim_and_actual_compiler(self):
        wrapper = self.folder / "rustc"
        wrapper.write_bytes(b"multicall compiler shim fixture")
        actual = self.folder / "toolchain/bin/rustc"
        actual.parent.mkdir(parents=True)
        actual.write_bytes(b"actual selected Rust compiler fixture")
        with mock.patch.object(DRIVER.shutil, "which", return_value=str(wrapper)), \
                mock.patch.object(DRIVER.subprocess, "check_output", side_effect=[
                    "rustc 1.90.0 (fixture)\n", str(actual.parent.parent) + "\n"]):
            compiler = DRIVER.compiler_identity("rustc")
        self.assertEqual(compiler["path"], str(wrapper))
        self.assertEqual(compiler["executable"], str(actual))
        self.assertEqual(compiler["sha256"], DRIVER.digest(actual))
        self.assertEqual(compiler["invocation_sha256"], DRIVER.digest(wrapper))

    def test_missing_compiler_acquires_only_private_pinned_prefix(self):
        prefix = self.folder / ("rust-" + DRIVER.ISOLATED_RUST_VERSION)
        compiler = dict(path=str(prefix / "bin/rustc"), version="rustc 1.90.0 (fixture)", sha256="f" * 64)
        with mock.patch.object(DRIVER, "compiler_identity", side_effect=[None, None, compiler]), \
                mock.patch.object(DRIVER.shutil, "which", return_value="/existing/bin/conda"), \
                mock.patch.object(DRIVER.subprocess, "run") as run:
            self.assertEqual(DRIVER.installation_compiler(self.folder), compiler)
        self.assertEqual(run.call_count, 1)
        self.assertEqual(run.call_args.args[0],
                         ["/existing/bin/conda", "create", "--yes", "--prefix", str(prefix),
                          "--override-channels", "--channel", "conda-forge", "rust=1.90.0"])
        self.assertTrue(run.call_args.kwargs["check"])
        self.assertEqual(run.call_args.kwargs["env"]["CONDA_PKGS_DIRS"], str(self.folder / ".conda-pkgs"))

    def test_incompatible_explicit_compiler_fails_without_acquisition(self):
        with mock.patch.object(DRIVER, "compiler_identity", return_value=None), \
                mock.patch.object(DRIVER.subprocess, "run") as run:
            with self.assertRaisesRegex(ValueError, "explicit --rustc"):
                DRIVER.installation_compiler(self.folder, rustc="/explicit/rustc")
        run.assert_not_called()

    def test_micromamba_on_path_is_supported_without_exported_manager(self):
        prefix = self.folder / ("rust-" + DRIVER.ISOLATED_RUST_VERSION)
        compiler = dict(path=str(prefix / "bin/rustc"), version="rustc 1.90.0 (fixture)", sha256="f" * 64)
        with mock.patch.dict(os.environ, {"MAMBA_EXE": "", "CONDA_EXE": ""}), \
                mock.patch.object(DRIVER, "compiler_identity", side_effect=[None, None, compiler]), \
                mock.patch.object(DRIVER.shutil, "which", side_effect=lambda name: "/existing/bin/micromamba" if name == "micromamba" else None), \
                mock.patch.object(DRIVER.subprocess, "run") as run:
            self.assertEqual(DRIVER.installation_compiler(self.folder), compiler)
        self.assertEqual(run.call_args.args[0][0], "/existing/bin/micromamba")

    def test_private_prefix_cannot_reuse_a_different_compiler_version(self):
        compiler = dict(path="/private/rustc", version="rustc 1.85.0 (fixture)", sha256="f" * 64)
        with mock.patch.object(DRIVER, "compiler_identity", side_effect=[None, compiler]), \
                mock.patch.object(DRIVER.subprocess, "run") as run:
            with self.assertRaisesRegex(ValueError, "does not match pinned"):
                DRIVER.installation_compiler(self.folder)
        run.assert_not_called()

    def test_cached_converter_reuse_does_not_consult_any_compiler(self):
        source = ROOT / "scripts/spliceai2_converter.rs"
        source_sha = DRIVER.digest(source)
        destination = self.folder / source_sha
        destination.mkdir()
        binary = destination / "spliceai2-converter"
        binary.write_bytes(b"completed converter fixture")
        binary.chmod(0o755)
        receipt = destination / "build.ready.json"
        receipt.write_text(json.dumps(dict(schema_version=1, source_sha256=source_sha,
                                          binary_sha256=DRIVER.digest(binary))))
        with mock.patch.object(DRIVER, "installation_compiler", side_effect=AssertionError("compiler-free reuse")):
            self.assertEqual(DRIVER.build_converter(self.folder), binary)
            binary.write_bytes(b"tampered converter fixture")
            with self.assertRaisesRegex(ValueError, "identity mismatch"):
                DRIVER.build_converter(self.folder)


class ProjectionTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.temporary = tempfile.TemporaryDirectory(prefix="spliceai2_projection_",
                                                    dir=os.environ.get("PRIVA_TEST_TMPDIR", ROOT / "tests"))
        cls.addClassCleanup(cls.temporary.cleanup)
        cls.folder = Path(cls.temporary.name)
        selected = os.environ.get("PRIVA_SPLICEAI2_INSTALL_BIN")
        if selected:
            cls.converter = Path(selected)
        else:
            if DRIVER.compiler_identity("rustc") is None:
                raise unittest.SkipTest("Use PRIVA_SPLICEAI2_INSTALL_BIN or a compatible rustc; unit tests do not acquire compilers")
            cls.converter = DRIVER.build_converter(cls.folder / "tools")
        cls.reference = os.environ.get("PRIVA_SPLICEAI2_REFERENCE_BIN")
        cls.genes = cls.folder / "genes.tsv"
        cls.genes.write_text(
            "gene_id\tchrom\tgene_start\tgene_end\tstrand\tgene_name\n"
            "PLUS_Z\tchr1\t10\t30\t+\tZETA\n"
            "PLUS_A\tchr1\t10\t30\t+\tALPHA\n"
            "PLUS_DUP\tchr1\t10\t30\t+\tALPHA\n"
            "MINUS\tchr1\t10\t30\t-\tMINUS\n")

    def invoke(self, rows, *, kind="snv", reference=False, header=False):
        command = [str(self.reference if reference else self.converter),
                   "spliceai2-vcf" if reference else "project",
                   "--genes", str(self.genes), "--chrom", "chr1", "--kind", kind,
                   "--input", "-", "--output", "-"]
        if not header:
            command.append("--no-header")
        content = "\t".join(native_columns()) + "\n" + "\n".join(rows) + "\n"
        return subprocess.run(command, input=content.encode(), capture_output=True, timeout=20)

    def test_source_order_gene_order_strand_and_score_lexemes(self):
        result = self.invoke([row(strand="-"), row(), row(pos="12")])
        self.assertEqual(result.returncode, 0, result.stderr)
        lines = result.stdout.decode().splitlines()
        self.assertEqual([line.split("|")[1] for line in lines], ["MINUS", "ALPHA", "ZETA", "ALPHA", "ZETA"])
        self.assertTrue(all("|0.2500|0.00|-0.0|0.00|-3|0|0|0" in line for line in lines))
        self.assertEqual(json.loads(result.stderr), dict(module="spliceai2-vcf", chrom="chr1", kind="snv",
                                                        input_rows=3, output_records=5))

    def test_intact_mnv_and_indel_alleles(self):
        result = self.invoke([row(ref="GG", alt="AA"), row(pos="12", ref="AC", alt="A", strand="-")], kind="indel")
        self.assertEqual(result.returncode, 0, result.stderr)
        lines = result.stdout.decode().splitlines()
        self.assertEqual([line.split("\t")[3:5] for line in lines], [["GG", "AA"], ["GG", "AA"], ["AC", "A"]])
        self.assertIn("SpliceAI=AA|ALPHA|", lines[0])

    def test_projection_header_preserves_official_v2_plugin_version_token(self):
        result = self.invoke([row()], header=True)
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertIn(b"SpliceAIv2.0", result.stdout)
        self.assertIn(b"annotation_only_no_ACMG_or_CADD_activation", result.stdout)

    def test_unprojected_missingness_duplicates_and_bad_scores_are_rejected(self):
        cases = [
            [row(jxn_loss_alt_score_2="")],
            [row(jxn_loss_alt_score_2="NaN")],
            [row(acceptor_gain_dist_0="2147483648")],
            [row(spliceai2_summary_score="0.10")],
            [row(), row()],
            [row(pos="12"), row(pos="11")],
            [row(ref="GG", alt="AA")],
        ]
        for rows in cases:
            with self.subTest(rows=rows):
                result = self.invoke(rows)
                self.assertEqual(result.returncode, 2)
                self.assertIn(b"native TSV line", result.stderr)

    def test_valid_streams_and_diagnostics_match_frozen_production_converter(self):
        if not self.reference:
            self.skipTest("Set PRIVA_SPLICEAI2_REFERENCE_BIN to the existing pinned production converter")
        for kind, rows in [
            ("snv", [row(strand="-"), row(), row(pos="12", acceptor_gain_delta_score_0="-0.25",
                                               donor_gain_delta_score_0="-0.1", donor_loss_delta_score_0="-0.2",
                                               acceptor_loss_delta_score_0="-0.3", spliceai2_summary_score="-0.1")]),
            ("indel", [row(ref="GG", alt="AA"), row(pos="12", ref="A", alt="ATT", strand="-")]),
        ]:
            with self.subTest(kind=kind):
                actual, expected = self.invoke(rows, kind=kind), self.invoke(rows, kind=kind, reference=True)
                self.assertEqual((actual.returncode, actual.stdout, actual.stderr),
                                 (expected.returncode, expected.stdout, expected.stderr))

    def test_failure_behavior_matches_frozen_production_converter(self):
        if not self.reference:
            self.skipTest("Set PRIVA_SPLICEAI2_REFERENCE_BIN to the existing pinned production converter")
        for rows in ([row(jxn_gain_ref_score_2="")], [row(acceptor_gain_delta_score_2="NaN")],
                     [row(), row()], [row(pos="12"), row(pos="11")]):
            with self.subTest(rows=rows):
                actual, expected = self.invoke(rows), self.invoke(rows, reference=True)
                self.assertEqual(actual.returncode, expected.returncode)
                self.assertEqual(actual.stdout, expected.stdout)
                # Executable name is the only intended error-diagnostic difference.
                self.assertEqual(actual.stderr.split(b": ", 1)[1], expected.stderr.split(b": ", 1)[1])


if __name__ == "__main__":
    unittest.main()

