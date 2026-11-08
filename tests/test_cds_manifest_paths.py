"""Nextflow CDS cache paths in the resource inventory; no scientific algorithm stubs.

The private task copies of the FASTA and its cache are tested in Rust
(rust/src/nextflow_task/tests.rs, cds_cache_copies_follow_the_detector_names).
"""
import gzip
import os
from pathlib import Path
import sys
import tempfile
import unittest


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "tests"))
import priva_task_cli as adapter
sys.path.insert(0, str(ROOT / "scripts"))
from alternative_start_codon import AlternativeStartCodonDetector


FASTA = b">ENST00000000233.1 cds chromosome:GRCh37:1:1:18:1\nATGAAAATGCCCATGTAA\n"


class CDSManifestPaths(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory(
            prefix="priva_cds_manifest_", dir=os.environ.get("PRIVA_TEST_TMPDIR", os.environ.get("TMPDIR")))
        self.addCleanup(self.temporary.cleanup)
        self.folder = Path(self.temporary.name).resolve()
        self.source = self.folder / "source"
        resource = self.folder / "resource.txt"
        resource.write_text("fixture\n")
        priva = self.folder / "priva"
        priva.write_text("#!/bin/sh\necho 'priva 0.1.0 (fixture)'\n")
        priva.chmod(0o755)
        self.config = {key: str(resource) for key in adapter.PRIORITY_KEYS}
        self.config.update(assembly="hg19", clinvar_gene_stat="baseline.pkl.gz",
                           repeat_region_file_name="repeat.bed.gz", dispensable_gene_list=str(resource),
                           priva_binary=str(priva))
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

    def test_runtime_and_adapter_agree_for_standard_plain_generic_and_parent_tokens(self):
        for index, relative in enumerate(("standard.fa.gz", "plain.fa", "generic.other.gz",
                                          "refs.fa.gz/selected.fa", "refs.fa.gz/selected.other.gz",
                                          "refs.fa.gz/standard.fa.gz")):
            with self.subTest(relative=relative):
                selected = self.fasta(relative)
                before = selected.read_bytes()
                detector = AlternativeStartCodonDetector(str(selected))
                cache = Path(detector.cache_pkl_path)
                self.config["cds_fasta_file"] = str(selected)
                resources = adapter.runtime_resources(self.config, self.source)
                self.assertEqual(resources["cds_parsed_cache"], cache)
                self.assertEqual(selected.read_bytes(), before)

    def test_old_whole_path_sidecar_never_overwrites_private_fasta_or_enters_inventory(self):
        for index, relative in enumerate(("refs.fa.gz/selected.fa", "refs.fa.gz/selected.other.gz")):
            with self.subTest(relative=relative):
                selected = self.fasta(relative)
                before = selected.read_bytes()
                old_cache = Path(str(selected).replace(".fa.gz", ".cds_cache.pkl.gz"))
                old_cache.parent.mkdir(parents=True, exist_ok=True)
                old_cache.write_bytes(b"old misplaced pickle fixture must not replace FASTA\n")
                self.config["cds_fasta_file"] = str(selected)
                resources = adapter.runtime_resources(self.config, self.source)
                self.assertNotIn("cds_parsed_cache", resources)
                self.assertEqual(selected.read_bytes(), before)
                self.assertTrue(old_cache.exists())

    def test_existing_nextflow_whitespace_guard_remains_active(self):
        selected = self.fasta("selected input's.fa")
        self.config["cds_fasta_file"] = str(selected)
        with self.assertRaises(adapter.ContractError):
            adapter.runtime_resources(self.config, self.source)


if __name__ == "__main__":
    unittest.main()
