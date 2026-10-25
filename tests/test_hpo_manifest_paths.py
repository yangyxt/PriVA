"""Manifest records the HPO file selected by the compatibility shell."""
import gzip
import hashlib
import importlib.util
import os
from pathlib import Path
import tempfile
import unittest


ROOT = Path(__file__).resolve().parents[1]
SPEC = importlib.util.spec_from_file_location("priva_task_hpo_paths", ROOT / "bin/priva_task.py")
adapter = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(adapter)


class HPOManifestPaths(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory(
            prefix="priva_hpo_manifest_", dir=os.environ.get("PRIVA_TEST_TMPDIR", os.environ.get("TMPDIR")))
        self.addCleanup(self.temporary.cleanup)
        self.folder = Path(self.temporary.name)
        self.source = self.folder / "source"
        self.fixture = self.folder / "resource.txt"
        self.fixture.write_text("routing fixture\n")
        self.config = {key: str(self.fixture) for key in adapter.PRIORITY_KEYS}
        self.config.update(assembly="hg19", clinvar_gene_stat="baseline.pkl.gz",
                           repeat_region_file_name="repeat.bed.gz", dispensable_gene_list=str(self.fixture))
        for relative in ("data/hpo/genes_to_phenotype.assertions.tsv.gz",
                         "data/alphamissense/alphamissense_mean_score.tsv", "data/ClinVar/baseline.pkl.gz",
                         "data/repeats/repeat.bed.gz", "data/hgnc/non_alt_loci_set.tsv"):
            path = self.source / relative
            path.parent.mkdir(parents=True, exist_ok=True)
            path.write_text("legacy routing fixture\n")
        self.default = self.source / "data/hpo/genes_to_phenotype.assertions.tsv.gz"
        self.selected = self.folder / "data_side_assertions.tsv.gz"
        with gzip.open(self.selected, "wt") as stream:
            stream.write("gene_symbol\tdisease_id\thpo_id\tfrequency\tevidence\treference\n"
                         "GENE\tOMIM:1\tHP:0000006\t-\tTAS\tPMID:1\n")

    def paths(self):
        return adapter.runtime_resources(self.config, self.source, "preannotated")

    def test_explicit_selected_hpo_identity_is_inventoried_without_legacy_substitution(self):
        self.config["hpo_assertions"] = str(self.selected)
        before = self.default.read_bytes()
        paths = self.paths()
        self.assertEqual(paths["hpo_assertions_used_by_shell"], self.selected)
        record = adapter.inventory(paths, checksum=True)["hpo_assertions_used_by_shell"]
        self.assertEqual(record["path"], str(self.selected))
        self.assertEqual(record["resolved_path"], str(self.selected.resolve()))
        self.assertEqual(record["sha256"], hashlib.sha256(self.selected.read_bytes()).hexdigest())
        self.assertEqual(self.default.read_bytes(), before)

    def test_missing_key_empty_and_null_keep_legacy_path(self):
        self.assertEqual(self.paths()["hpo_assertions_used_by_shell"], self.default)
        for value in (None, "", "null"):
            with self.subTest(value=value):
                self.config["hpo_assertions"] = value
                self.assertEqual(self.paths()["hpo_assertions_used_by_shell"], self.default)

    def test_explicit_missing_hpo_fails_with_existing_legacy_file(self):
        self.config["hpo_assertions"] = str(self.folder / "missing.tsv.gz")
        with self.assertRaisesRegex(adapter.ContractError, "hpo_assertions_used_by_shell: missing"):
            self.paths()
        self.assertTrue(self.default.exists())

    def test_selected_file_symlink_records_the_actual_target_identity(self):
        selected_link = self.folder / "selected_link.tsv.gz"
        selected_link.symlink_to(self.selected)
        self.config["hpo_assertions"] = str(selected_link)
        record = adapter.inventory(self.paths(), checksum=True)["hpo_assertions_used_by_shell"]
        self.assertEqual(record["path"], str(selected_link))
        self.assertEqual(record["resolved_path"], str(self.selected.resolve()))
        self.assertEqual(record["sha256"], hashlib.sha256(self.selected.read_bytes()).hexdigest())

    def test_existing_nextflow_path_safety_contract_is_preserved(self):
        invalid = self.folder / "quoted HPO's path.tsv.gz"
        invalid.write_bytes(self.selected.read_bytes())
        for value in (str(invalid), False, True, "relative.tsv.gz"):
            with self.subTest(value=value):
                self.config["hpo_assertions"] = value
                with self.assertRaises(adapter.ContractError):
                    self.paths()


if __name__ == "__main__":
    unittest.main()
