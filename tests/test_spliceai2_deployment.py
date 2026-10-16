"""Synthetic deployment tests; real-model, full-cache and VEP acceptance are separate."""
import argparse
import contextlib
import gzip
import importlib.util
import io
import json
import os
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest
from unittest import mock

ROOT = Path(__file__).resolve().parents[1]
SPEC = importlib.util.spec_from_file_location("spliceai2_deployment", ROOT / "scripts/spliceai2_to_vep_vcf.py")
DEPLOY = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(DEPLOY)


def native_columns():
    columns = ["chrom", "pos", "ref", "alt", "strand"]
    for effect in ("donor_gain", "donor_loss", "acceptor_gain", "acceptor_loss", "jxn_gain", "jxn_loss"):
        fields = ["delta_score", "ref_score", "alt_score"] + (["donor_dist", "acceptor_dist"] if effect.startswith("jxn") else ["dist"])
        columns.extend(f"{effect}_{field}_{index}" for field in fields for index in range(3))
    return columns + ["spliceai2_summary_score"]


class GeneFixture:
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory(prefix="spliceai2_deploy_", dir=ROOT / "tests")
        self.addCleanup(self.temp.cleanup)
        self.folder = Path(self.temp.name)
        self.regions = self.folder / "regions.tsv"
        self.regions.write_text("gene_id\tchrom\tgene_start\tgene_end\tstrand\nGENE1.1\tchr1\t10\t30\t+\nGENE2.2\tchr1\t10\t30\t-\n")
        self.gtf = self.folder / "gencode.v48.annotation.gtf.gz"
        content = '##description: GRCh38, version 48\n'
        for gene, name, strand in (("GENE1.1", "FIXTURE_A", "+"), ("GENE2.2", "FIXTURE_B", "-")):
            content += f'chr1\tfixture\tgene\t10\t30\t.\t{strand}\t.\tgene_id "{gene}"; gene_name "{name}"; gene_type "protein_coding"; tag "one"; tag "two";\n'
        self.gtf.write_bytes(gzip.compress(content.encode()))
        self.context = self.folder / "context.tsv"

    def build(self):
        return DEPLOY.gene_context(self.regions, self.gtf, self.context)


class GeneContextTests(GeneFixture, unittest.TestCase):
    def test_exact_gtf_context_preserves_both_model_strands(self):
        report = self.build()
        self.assertEqual(report["source_coordinates"], "1-based-closed")
        self.assertEqual(report["gene_id_matching"], "exact")
        self.assertEqual(self.context.read_text().splitlines()[1:], [
            "GENE1.1\tchr1\t10\t30\t+\tFIXTURE_A", "GENE2.2\tchr1\t10\t30\t-\tFIXTURE_B"])

    def test_zero_based_only_when_every_bound_matches_gtf(self):
        self.regions.write_text(self.regions.read_text().replace("\t10\t30", "\t9\t30"))
        self.assertEqual(self.build()["source_coordinates"], "0-based-half-open")
        self.assertIn("\t10\t30\t", self.context.read_text())

    def test_mixed_coordinates_rejected(self):
        self.regions.write_text(self.regions.read_text().replace("\t10\t30", "\t9\t30", 1))
        with self.assertRaisesRegex(ValueError, "mixed-coordinate"):
            self.build()

    def test_no_silent_gene_version_stripping(self):
        self.regions.write_text(self.regions.read_text().replace("GENE1.1", "GENE1"))
        with self.assertRaisesRegex(ValueError, "No exact or unique GENCODE48 stable gene ID"):
            self.build()

    def test_gene_bounds_and_strand_mismatches_rejected(self):
        original = self.regions.read_text()
        for changed in (original.replace("\t10\t30", "\t8\t30", 1), original.replace("\t+\n", "\t-\n", 1)):
            self.regions.write_text(changed)
            with self.assertRaises(ValueError):
                self.build()
            self.assertFalse(self.context.exists())

    def test_wrong_gtf_release_rejected(self):
        self.gtf.write_bytes(gzip.compress(gzip.decompress(self.gtf.read_bytes()).replace(b"version 48", b"version 50")))
        with self.assertRaisesRegex(ValueError, "GENCODE version 48"):
            self.build()

    def test_duplicate_gene_id_rejected(self):
        self.regions.write_text(self.regions.read_text() + "GENE1.1\tchr1\t10\t30\t+\n")
        with self.assertRaisesRegex(ValueError, "Duplicate scored-region gene"):
            self.build()

    def test_same_symbol_opposite_strand_gene_ids_on_one_chromosome_rejected(self):
        self.gtf.write_bytes(gzip.compress(gzip.decompress(self.gtf.read_bytes()).replace(b"FIXTURE_B", b"FIXTURE_A")))
        with self.assertRaisesRegex(ValueError, "Opposite-strand genes share a symbol"):
            self.build()

    def test_same_symbol_same_strand_retains_every_gene_context(self):
        self.regions.write_text(self.regions.read_text().replace("\t-\n", "\t+\n"))
        self.gtf.write_bytes(gzip.compress(gzip.decompress(self.gtf.read_bytes()).replace(b"FIXTURE_B", b"FIXTURE_A").replace(b"\t-\t", b"\t+\t")))
        report = self.build()
        self.assertEqual(report["same_strand_symbol_groups"], 1)
        self.assertEqual(len(self.context.read_text().splitlines()), 3)

    def test_publisher_stable_ids_resolve_unique_versioned_ids_without_losing_par_y(self):
        self.regions.write_text(self.regions.read_text().replace("GENE1.1", "ENSG00000186092").replace("GENE2.2", "ENSG00000124333_PAR_Y"))
        self.gtf.write_bytes(gzip.compress(gzip.decompress(self.gtf.read_bytes()).replace(b"GENE1.1", b"ENSG00000186092.7").replace(b"GENE2.2", b"ENSG00000124333.9_PAR_Y")))
        report = self.build()
        self.assertEqual(report["stable_id_matches"], 2)
        self.assertEqual(report["gene_id_matching"], "unique-GENCODE48-stable-ID-with-exact-context")
        self.assertIn("ENSG00000186092.7\t", self.context.read_text())
        self.assertIn("ENSG00000124333.9_PAR_Y\t", self.context.read_text())

    def test_ambiguous_stable_id_versions_rejected(self):
        self.regions.write_text(self.regions.read_text().replace("GENE1.1", "ENSG00000186092"))
        content = gzip.decompress(self.gtf.read_bytes()).decode().replace("GENE1.1", "ENSG00000186092.7")
        content += content.splitlines()[1].replace("ENSG00000186092.7", "ENSG00000186092.8") + "\n"
        self.gtf.write_bytes(gzip.compress(content.encode()))
        with self.assertRaisesRegex(ValueError, "aliases/ambiguous versions"):
            self.build()

    def test_partial_symlink_target_rejected(self):
        partial = self.folder / "source.partial"
        partial.write_text("unfinished")
        alias = self.folder / "complete.tsv"
        alias.symlink_to(partial)
        with self.assertRaisesRegex(ValueError, "Incomplete input target"):
            DEPLOY.complete_file(alias)


class DeploymentTests(GeneFixture, unittest.TestCase):
    def setUp(self):
        super().setUp()
        # The driver fixture pins a synthetic plugin identity. Actual official
        # plugin parsing/allele matching is exercised by the separate Perl test.
        self.plugin = self.folder / "SpliceAI.pm"
        self.plugin.write_text("synthetic identity fixture; not an executed plugin\n")
        self.addCleanup(mock.patch.stopall)
        mock.patch.object(DEPLOY, "PLUGIN_SHA256", DEPLOY.digest(self.plugin)).start()
        self.manifest = self.folder / "confirmed.json"
        self.files = []
        for kind in ("snv", "indel"):
            path = self.folder / f"SpliceAI2_GRCh38_{kind}_chr1.tsv.gz"
            columns = native_columns()
            rows = []
            for position, score in ((11, "0.25"), (12, "0.09")):
                values = {name: ("0" if "dist_" in name else "0.00") for name in columns}
                values.update(chrom="chr1", pos=str(position), ref="A" if kind == "snv" else "AC",
                              alt="G" if kind == "snv" else "A", strand="+",
                              acceptor_gain_delta_score_0=score, donor_gain_delta_score_0="-0.0",
                              acceptor_gain_dist_0="-3", spliceai2_summary_score=score)
                rows.append("\t".join(values[name] for name in columns))
            path.write_bytes(gzip.compress(("\t".join(columns) + "\n" + "\n".join(rows) + "\n").encode()))
            self.files.append({"kind": kind, "chrom": "chr1", **self.item(path)})
        self.document = {"schema_version": 1, "status": "downloads_complete", "confirmed_by": "synthetic fixture",
                         "resource": "SpliceAI2", "version": "2.0", "assembly": "GRCh38", "source_revision": DEPLOY.REVISION,
                         "files": self.files, "scored_regions": self.item(self.regions),
                         "gencode_gtf": {"release": "GENCODE48", **self.item(self.gtf)},
                         "vep_spliceai_plugin": self.item(self.plugin)}
        self.manifest.write_text(json.dumps(self.document))
        self.output = self.folder / "prepared"
        self.args = argparse.Namespace(manifest=self.manifest, assembly="GRCh38", output_dir=self.output,
                                       receipt=self.output / "adapter.ready.json", priva_bin=os.environ.get("PRIVA_SPLICEAI2_BIN", sys.executable),
                                       threads=1, validation_subset=True)

    @staticmethod
    def item(path):
        return {"path": str(path), "bytes": path.stat().st_size, "sha256": DEPLOY.digest(path)}

    def prepare(self):
        with contextlib.redirect_stdout(io.StringIO()):
            return DEPLOY.prepare(self.args)

    def require_converter(self):
        if not os.environ.get("PRIVA_SPLICEAI2_BIN"):
            self.skipTest("Set PRIVA_SPLICEAI2_BIN to the built real Rust priva for deployment integration")

    def test_unconfirmed_input_rejected_before_sources(self):
        self.document["status"] = "downloading"
        self.document["files"][0]["path"] = str(self.folder / "does-not-exist")
        self.manifest.write_text(json.dumps(self.document))
        with self.assertRaisesRegex(ValueError, "confirmed, complete"):
            self.prepare()
        self.assertFalse(self.output.exists())

    def test_hg19_route_rejected(self):
        self.args.assembly = "hg19"
        with self.assertRaisesRegex(ValueError, "hg19 remains v1"):
            self.prepare()

    def test_boolean_manifest_version_rejected(self):
        self.document["schema_version"] = True
        self.manifest.write_text(json.dumps(self.document))
        with self.assertRaisesRegex(ValueError, "confirmed, complete"):
            self.prepare()

    def test_partial_or_symlink_output_directory_rejected(self):
        self.output.mkdir()
        alias = self.folder / "output-alias"
        alias.symlink_to(self.output, target_is_directory=True)
        for path, message in ((alias, "symlinks"), (self.folder / "output.PARTIAL", "Incomplete/noncanonical")):
            self.args.output_dir = path
            self.args.receipt = path / "adapter.ready.json"
            with self.assertRaisesRegex(ValueError, message):
                self.prepare()
        self.assertFalse((self.output / "adapter.ready.json").exists())

    def test_subset_cannot_claim_complete_published_coverage(self):
        self.args.validation_subset = False
        with self.assertRaisesRegex(ValueError, "all 47 published shards"):
            self.prepare()

    def test_source_checksum_mismatch_rejected(self):
        self.document["files"][0]["sha256"] = "0" * 64
        self.manifest.write_text(json.dumps(self.document))
        with self.assertRaisesRegex(ValueError, "SHA256 mismatch"):
            self.prepare()

    def test_real_rust_bgzip_tabix_subset_retains_low_scores_and_originals(self):
        self.require_converter()
        originals = {Path(item["path"]): Path(item["path"]).read_bytes() for item in self.files}
        report = self.prepare()
        self.assertEqual(report["status"], "validated_subset_annotation_only")
        self.assertFalse(report["scientific_qualified"])
        self.assertEqual(report["thresholds"], {"high_recall": .1, "balanced": .25, "high_precision": .5})
        for kind in ("snv", "indel"):
            artifact = Path(report["outputs"][kind]["path"])
            rows = subprocess.check_output(["tabix", str(artifact), "chr1:10-20"], text=True).splitlines()
            self.assertEqual(len(rows), 2)
            self.assertIn("|FIXTURE_A|0.25|0.00|-0.0|0.00|-3|0|0|0", rows[0])
            self.assertIn("|FIXTURE_A|0.09|", rows[1])
            self.assertTrue(all("FIXTURE_B" not in row for row in rows))
            self.assertTrue(all(row.split("\t")[7].startswith("SpliceAI=") and ";" not in row.split("\t")[7] for row in rows))
        for path, content in originals.items():
            self.assertEqual(path.read_bytes(), content)
        self.assertEqual(list(self.output.glob(".building-*")), [])

    def test_changed_source_during_projection_never_publishes_receipt(self):
        self.require_converter()
        actual = DEPLOY.convert_kind
        target = Path(self.files[0]["path"])
        def mutate(*args, **kwargs):
            result = actual(*args, **kwargs)
            stat = target.stat()
            os.utime(target, ns=(stat.st_atime_ns, stat.st_mtime_ns + 1_000_000_000))
            return result
        with mock.patch.object(DEPLOY, "convert_kind", side_effect=mutate):
            with self.assertRaisesRegex(ValueError, "Input changed while preparing"):
                self.prepare()
        self.assertFalse(self.args.receipt.exists())
        self.assertEqual(list(self.output.glob("*.vcf.gz*")), [])

    def test_real_projection_coalesces_identical_symbol_context_without_losing_predictions(self):
        self.require_converter()
        self.regions.write_text(self.regions.read_text().replace("\t-\n", "\t+\n"))
        self.gtf.write_bytes(gzip.compress(gzip.decompress(self.gtf.read_bytes()).replace(b"FIXTURE_B", b"FIXTURE_A").replace(b"\t-\t", b"\t+\t")))
        self.document.update(scored_regions=self.item(self.regions),
                             gencode_gtf={"release": "GENCODE48", **self.item(self.gtf)})
        self.manifest.write_text(json.dumps(self.document))
        report = self.prepare()
        self.assertEqual(report["gene_context"]["genes"], 2)
        self.assertEqual(report["gene_context"]["same_strand_symbol_groups"], 1)
        for kind in ("snv", "indel"):
            rows = subprocess.check_output(["tabix", report["outputs"][kind]["path"], "chr1:10-20"], text=True).splitlines()
            self.assertEqual(len(rows), 2)
            self.assertEqual(report["statistics"][kind][0]["input_rows"], 2)
            self.assertEqual(report["statistics"][kind][0]["output_records"], 2)
            self.assertIn("|FIXTURE_A|0.25|", rows[0])
            self.assertIn("|FIXTURE_A|0.09|", rows[1])

    def test_full_47_shard_stream_has_one_header_and_exact_indexed_contigs(self):
        self.require_converter()
        original_regions = self.regions.read_text().splitlines()
        original_gtf = gzip.decompress(self.gtf.read_bytes()).decode().splitlines()
        region_lines, gtf_lines = [original_regions[0]], [original_gtf[0]]
        templates = {item["kind"]: gzip.decompress(Path(item["path"]).read_bytes()).decode() for item in self.files}
        shards = []
        for chrom in DEPLOY.CHROMS:
            for line in original_regions[1:]:
                region_lines.append(line.replace("chr1", chrom).replace("GENE", f"{chrom}_GENE"))
            for line in original_gtf[1:]:
                gtf_lines.append(line.replace("chr1", chrom).replace("GENE", f"{chrom}_GENE"))
            for kind in ("snv", "indel"):
                if (kind, chrom) not in DEPLOY.PUBLISHED:
                    continue
                path = self.folder / f"SpliceAI2_GRCh38_{kind}_{chrom}.tsv.gz"
                path.write_bytes(gzip.compress(templates[kind].replace("chr1", chrom).encode()))
                shards.append({"kind": kind, "chrom": chrom, **self.item(path)})
        self.regions.write_text("\n".join(region_lines) + "\n")
        self.gtf.write_bytes(gzip.compress(("\n".join(gtf_lines) + "\n").encode()))
        self.document.update(files=shards, scored_regions=self.item(self.regions),
                             gencode_gtf={"release": "GENCODE48", **self.item(self.gtf)})
        self.manifest.write_text(json.dumps(self.document))
        self.args.validation_subset = False
        report = self.prepare()
        self.assertEqual(report["status"], "prepared_annotation_only")
        self.assertEqual(report["coverage"], "complete_published")
        self.assertFalse(report["scientific_qualified"])
        self.assertEqual(len(report["source_files"]), 47)
        for kind in ("snv", "indel"):
            expected = DEPLOY.CHROMS if kind == "snv" else DEPLOY.CHROMS[:-1]
            artifact = Path(report["outputs"][kind]["path"])
            self.assertEqual(subprocess.check_output(["tabix", "-l", str(artifact)], text=True).splitlines(), expected)
            with gzip.open(artifact, "rt") as stream:
                lines = stream.readlines()
            self.assertEqual(sum(line.startswith("#CHROM") for line in lines), 1)
            self.assertEqual(sum(not line.startswith("#") for line in lines), len(expected) * 2)
            for chrom in expected:
                rows = subprocess.check_output(["tabix", str(artifact), f"{chrom}:10-20"], text=True).splitlines()
                self.assertEqual(len(rows), 2)

    def test_corrupt_gzip_crc_never_publishes_output_or_receipt(self):
        self.require_converter()
        raw = Path(self.files[0]["path"])
        data = bytearray(raw.read_bytes())
        data[-8] ^= 1  # Intact data rows, corrupted gzip CRC trailer.
        raw.write_bytes(data)
        self.document["files"][0].update(self.item(raw))
        self.manifest.write_text(json.dumps(self.document))
        with self.assertRaisesRegex(ValueError, "Source decompression failed"):
            self.prepare()
        self.assertFalse(self.args.receipt.exists())
        self.assertEqual(list(self.output.glob("*.vcf.gz*")), [])

    def test_rust_failure_never_publishes_output_or_receipt(self):
        self.require_converter()
        raw = Path(self.files[0]["path"])
        raw.write_bytes(gzip.compress(gzip.decompress(raw.read_bytes()).replace(b"acceptor_gain_dist_0", b"unexpected_dist")))
        self.document["files"][0].update(self.item(raw))
        self.manifest.write_text(json.dumps(self.document))
        with self.assertRaisesRegex(ValueError, "Rust projection failed"):
            self.prepare()
        self.assertFalse(self.args.receipt.exists())
        self.assertEqual(list(self.output.glob("*.vcf.gz*")), [])

    def test_existing_prepared_files_never_overwritten(self):
        self.require_converter()
        self.prepare()
        before = self.args.receipt.read_bytes()
        with self.assertRaisesRegex(ValueError, "Refusing to replace"):
            self.prepare()
        self.assertEqual(self.args.receipt.read_bytes(), before)


if __name__ == "__main__":
    unittest.main()
