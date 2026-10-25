"""Nextflow CDS cache paths/inventory/private copies; no scientific algorithm stubs."""
import gzip
import importlib.util
import os
from pathlib import Path
import sys
import tempfile
import unittest


ROOT = Path(__file__).resolve().parents[1]
SPEC = importlib.util.spec_from_file_location("priva_task_cds_paths", ROOT / "bin/priva_task.py")
adapter = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(adapter)
sys.path.insert(0, str(ROOT / "scripts"))
from alternative_start_codon import AlternativeStartCodonDetector


FASTA = b">ENST00000000233.1 cds chromosome:GRCh37:1:1:18:1\nATGAAAATGCCCATGTAA\n"


class CDSManifestPaths(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory(
            prefix="priva_cds_manifest_", dir=os.environ.get("PRIVA_TEST_TMPDIR", os.environ.get("TMPDIR")))
        self.addCleanup(self.temporary.cleanup)
        self.folder = Path(self.temporary.name)
        self.source = self.folder / "source"
        resource = self.folder / "resource.txt"
        resource.write_text("fixture\n")
        self.config = {key: str(resource) for key in adapter.PRIORITY_KEYS}
        self.config.update(assembly="hg19", clinvar_gene_stat="baseline.pkl.gz",
                           repeat_region_file_name="repeat.bed.gz", dispensable_gene_list=str(resource))
        for relative in ("data/hpo/genes_to_phenotype.assertions.tsv.gz",
                         "data/alphamissense/alphamissense_mean_score.tsv", "data/ClinVar/baseline.pkl.gz",
                         "data/repeats/repeat.bed.gz", "data/hgnc/non_alt_loci_set.tsv"):
            path = self.source / relative
            path.parent.mkdir(parents=True, exist_ok=True)
            path.write_text("legacy curated source fixture\n")

    def fasta(self, relative):
        path = self.folder / relative
        path.parent.mkdir(parents=True, exist_ok=True)
        if path.name.endswith(".gz"):
            with gzip.open(path, "wb") as handle:
                handle.write(FASTA)
        else:
            path.write_bytes(FASTA)
        return path

    def private_config(self, selected, tag):
        bundle = self.folder / ("bundle_" + tag)
        bundle.mkdir()
        (bundle / "family.ped").write_text("F SAMPLE 0 0 1 2\n")
        (bundle / "legacy.yaml").write_text(adapter.yaml.safe_dump({"assembly": "hg19", "cds_fasta_file": str(selected)}))
        task = self.folder / ("task_" + tag)
        task.mkdir()
        previous = Path.cwd()
        try:
            os.chdir(task)
            path, _ = adapter.private_config(bundle, {"base": "public_fixture", "assembly": "hg19", "mode": "preannotated"}, "prioritization", 1)
            return adapter.yaml.safe_load(path.read_text())
        finally:
            os.chdir(previous)

    def test_runtime_and_adapter_agree_for_standard_plain_generic_and_parent_tokens(self):
        for index, relative in enumerate(("standard.fa.gz", "plain.fa", "generic.other.gz",
                                          "refs.fa.gz/selected.fa", "refs.fa.gz/selected.other.gz",
                                          "refs.fa.gz/standard.fa.gz")):
            with self.subTest(relative=relative):
                selected = self.fasta(relative)
                before = selected.read_bytes()
                detector = AlternativeStartCodonDetector(str(selected))
                cache = Path(detector.cache_pkl_path)
                self.assertEqual(adapter.cds_cache_path(selected), cache)
                self.config["cds_fasta_file"] = str(selected)
                resources = adapter.runtime_resources(self.config, self.source, "preannotated")
                self.assertEqual(resources["cds_parsed_cache"], cache)
                configured = self.private_config(selected, str(index))
                private_fasta = Path(configured["cds_fasta_file"])
                self.assertEqual(private_fasta.read_bytes(), before)
                self.assertEqual(adapter.cds_cache_path(private_fasta).read_bytes(), cache.read_bytes())
                self.assertEqual(selected.read_bytes(), before)
                private_reader = AlternativeStartCodonDetector(str(private_fasta))
                self.assertEqual(private_reader.find_alternative_start("ENST00000000233"),
                                 detector.find_alternative_start("ENST00000000233"))
                self.assertEqual(private_fasta.read_bytes(), before)

    def test_old_whole_path_sidecar_never_overwrites_private_fasta_or_enters_inventory(self):
        for index, relative in enumerate(("refs.fa.gz/selected.fa", "refs.fa.gz/selected.other.gz")):
            with self.subTest(relative=relative):
                selected = self.fasta(relative)
                before = selected.read_bytes()
                old_cache = Path(str(selected).replace(".fa.gz", ".cds_cache.pkl.gz"))
                old_cache.parent.mkdir(parents=True, exist_ok=True)
                old_cache.write_bytes(b"old misplaced pickle fixture must not replace FASTA\n")
                self.config["cds_fasta_file"] = str(selected)
                resources = adapter.runtime_resources(self.config, self.source, "preannotated")
                self.assertNotIn("cds_parsed_cache", resources)
                configured = self.private_config(selected, "old_" + str(index))
                private_fasta = Path(configured["cds_fasta_file"])
                self.assertEqual(private_fasta.read_bytes(), before)
                self.assertFalse(adapter.cds_cache_path(private_fasta).exists())
                self.assertEqual(selected.read_bytes(), before)
                self.assertTrue(old_cache.exists())

    def test_existing_nextflow_whitespace_guard_remains_active(self):
        selected = self.fasta("selected input's.fa")
        self.config["cds_fasta_file"] = str(selected)
        with self.assertRaises(adapter.ContractError):
            adapter.runtime_resources(self.config, self.source, "preannotated")


if __name__ == "__main__":
    unittest.main()
