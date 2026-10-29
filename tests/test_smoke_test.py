"""Boundary tests for the prepared-fixture smoke runner, not model acceptance."""
import argparse
import contextlib
import csv
import gzip
import importlib.util
import io
from pathlib import Path
import pickle
import tempfile
import unittest

import yaml

ROOT = Path(__file__).resolve().parents[1]
spec = importlib.util.spec_from_file_location("priva_smoke", ROOT / "example/smoke_test.py")
SMOKE = importlib.util.module_from_spec(spec)
spec.loader.exec_module(SMOKE)


class SmokeBoundaries(unittest.TestCase):
    def setUp(self):
        self.workspace = tempfile.TemporaryDirectory()
        self.addCleanup(self.workspace.cleanup)
        self.folder = Path(self.workspace.name)
        opaque = self.folder / "resource"
        opaque.write_bytes(b"Bundled cache existence/size check only.\n")
        self.config = {key: str(opaque) for key in SMOKE.REQUIRED_RESOURCES}
        self.vcf = self.folder / "prepared.vcf.gz"
        with gzip.open(self.vcf, "wt") as stream:
            stream.write('##fileformat=VCFv4.2\n##INFO=<ID=CSQ,Number=.,Type=String,Description="Format: Allele|Feature">\n')
            stream.write('#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\nchr1\t10\t.\tA\tC\t.\tPASS\tCSQ=C|ENST1\n')
        Path(str(self.vcf) + ".tbi").write_bytes(b"Existence check; real indexes exercised by end-to-end smoke runs.\n")
        self.cadd = self.folder / "prepared.cadd.tsv"
        self.cadd.write_text('#Chrom\tPos\tRef\tAlt\tFeatureID\tPHRED\nchr1\t10\tA\tC\tENST1\t12\n')
        self.config.update(smoke_fixture_only=True, assembly="hg19", ped_file="",
                           input_vcf=str(self.vcf), hub_cadd_file=str(self.cadd),
                           repeat_region_file_name="ucsc.hg19.fasta.repeat.bed.gz")
        self.source = self.folder / "case.yaml"

    def arguments(self, **kwargs):
        self.source.write_text(yaml.safe_dump(self.config))
        values = dict(assembly="hg19", config=self.source, output_dir=self.folder / "run", cores=2)
        values.update(kwargs)
        return argparse.Namespace(**values)

    def test_prepare_keeps_sources_unchanged_and_writes_only_external_paths(self):
        old = (self.vcf.read_bytes(), self.cadd.read_bytes())
        prepared = SMOKE.prepare(self.arguments())
        config = SMOKE.load_config(prepared)
        self.assertEqual(old, (self.vcf.read_bytes(), self.cadd.read_bytes()))
        self.assertEqual(Path(config["input_vcf"]).read_bytes(), old[0])
        self.assertEqual(Path(config["input_vcf"]).parent, prepared.parent / "input")
        self.assertEqual(Path(config["output_dir"]), prepared.parent / "output")
        self.assertNotIn(str(ROOT), config["output_dir"])

    def test_wrong_assembly_is_rejected_before_creating_output(self):
        args = self.arguments(assembly="hg38")
        with self.assertRaisesRegex(ValueError, "assembly"):
            SMOKE.prepare(args)
        self.assertFalse(args.output_dir.exists())

    def test_existing_output_is_preserved(self):
        args = self.arguments()
        args.output_dir.mkdir()
        marker = args.output_dir / "keep.txt"
        marker.write_text("Earlier smoke results must remain intact.")
        with self.assertRaisesRegex(ValueError, "not empty"):
            SMOKE.prepare(args)
        self.assertEqual(marker.read_text(), "Earlier smoke results must remain intact.")

    def test_checkout_output_is_rejected(self):
        with self.assertRaisesRegex(ValueError, "outside the source checkout"):
            SMOKE.prepare(self.arguments(output_dir=ROOT / "example/smoke_runs"))

    def test_missing_required_input_is_not_replaced(self):
        self.config["clinvar_aa_stat"] = str(self.folder / "absent.pkl")
        args = self.arguments()
        with self.assertRaisesRegex(ValueError, "Missing or empty"):
            SMOKE.prepare(args)
        self.assertFalse(args.output_dir.exists())

    def test_missing_cadd_allele_is_rejected_without_chromosome_rewriting(self):
        self.cadd.write_text(self.cadd.read_text().replace("chr1\t", "1\t"))
        args = self.arguments()
        with self.assertRaisesRegex(ValueError, "cover all input alleles"):
            SMOKE.prepare(args)
        self.assertFalse(args.output_dir.exists())

    def test_phred_only_table_is_rejected_without_feature_id(self):
        self.cadd.write_text('#Chrom\tPos\tRef\tAlt\tPHRED\nchr1\t10\tA\tC\t12\n')
        with self.assertRaisesRegex(ValueError, "FeatureID"):
            SMOKE.input_check(self.config)

    def test_incomplete_csq_is_rejected(self):
        with gzip.open(self.vcf, "wt") as stream:
            stream.write('##fileformat=VCFv4.2\n##INFO=<ID=CSQ,Number=.,Type=String,Description="Format: Allele|Feature">\n')
            stream.write('#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\nchr1\t10\t.\tA\tC\t.\tPASS\tCSQ=C\n')
        with self.assertRaisesRegex(ValueError, "width mismatch"):
            SMOKE.input_check(self.config)

    def test_production_config_cannot_select_smoke_workflow_accidentally(self):
        self.config.pop("smoke_fixture_only")
        self.source.write_text(yaml.safe_dump(self.config))
        with self.assertRaisesRegex(ValueError, "smoke_fixture_only"):
            SMOKE.load_config(self.source)

    def test_bundled_clinvar_aa_fixtures_load_without_main_module_factories(self):
        for assembly in ("hg19", "hg38"):
            path = ROOT / f"example/smoke_data/{assembly}/resources/clinvar.{assembly}.demo20.vep.aa_change.stats.pkl"
            with path.open("rb") as stream:
                cache = pickle.load(stream)
            self.assertIs(type(cache), dict)
            self.assertTrue(cache)

    def completed_outputs(self):
        path = SMOKE.prepare(self.arguments())
        config = SMOKE.load_config(path)
        base = Path(config["output_dir"]) / "demo20.hg19.anno.filtered"
        identity = {"chrom": "chr1", "pos": "10", "ref": "A", "alt": "C", "Feature": "ENST1", "BIOTYPE": "protein_coding"}
        other = dict(identity, Feature="ENST2", BIOTYPE="processed_transcript")
        for suffix, rows in [
            (".prepared.tsv", [identity, other]),
            (".tsv", [dict(identity, ACMG_criteria="PM2", ACMG_class="Uncertain_significance", ACMG_quant_score="0.5", sort_index="0.5", variant_rank="1")]),
            (".acmg.tsv", [{"PVS1": "0", "PM2": "1"}, {"PVS1": "0", "PM2": "1"}]),
        ]:
            target = base.parent / "demo20.hg19.anno.prepared.tsv" if suffix == ".prepared.tsv" else Path(str(base) + suffix)
            with target.open("w") as stream:
                writer = csv.DictWriter(stream, fieldnames=list(rows[0]), delimiter="\t")
                writer.writeheader()
                writer.writerows(rows)
        Path(str(base) + ".vcf.gz").write_bytes(self.vcf.read_bytes())
        Path(str(base) + ".cadd.tsv").write_bytes(self.cadd.read_bytes())
        return argparse.Namespace(config=path), base

    def test_matrix_keeps_noncoding_rows_while_ranking_keeps_protein_coding(self):
        args, _ = self.completed_outputs()
        output = io.StringIO()
        with contextlib.redirect_stdout(output):
            SMOKE.validate(args)
        self.assertIn("1 ranked and 2 ACMG rows", output.getvalue())

    def test_lost_or_changed_protein_transcript_is_reported(self):
        args, base = self.completed_outputs()
        table = Path(str(base) + ".tsv")
        table.write_text(table.read_text().replace("ENST1", "ENST_MISMATCH"))
        with self.assertRaisesRegex(ValueError, "FeatureID rows"):
            SMOKE.validate(args)

    def test_nonfinite_ranking_score_is_reported(self):
        args, base = self.completed_outputs()
        table = Path(str(base) + ".tsv")
        table.write_text(table.read_text().replace("0.5", "nan"))
        with self.assertRaisesRegex(ValueError, "Nonfinite score"):
            SMOKE.validate(args)


if __name__ == "__main__":
    unittest.main()
