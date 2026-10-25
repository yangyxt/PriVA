"""Inventory and readiness use the gene-statistics file selected by the shell."""
import contextlib
import gzip
import hashlib
import importlib.util
import io
import json
import os
from pathlib import Path
import pickle
import sys
import tempfile
import unittest
from unittest.mock import patch

import yaml


ROOT = Path(__file__).resolve().parents[1]


def load_module(name, relative):
    spec = importlib.util.spec_from_file_location(name, ROOT / relative)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


adapter = load_module("priva_task_clinvar_paths", "bin/priva_task.py")
readiness = load_module("readiness_clinvar_paths", "scripts/check_cluster_readiness.py")


@contextlib.contextmanager
def working_directory(path):
    previous = Path.cwd()
    try:
        os.chdir(path)
        yield
    finally:
        os.chdir(previous)


class ClinVarManifestPaths(unittest.TestCase):
    def setUp(self):
        temporary = tempfile.TemporaryDirectory(
            prefix="priva_clinvar_manifest_", dir=os.environ.get("PRIVA_TEST_TMPDIR", os.environ.get("TMPDIR")))
        self.addCleanup(temporary.cleanup)
        self.folder = Path(temporary.name)
        self.source = self.folder / "source"
        resource = self.folder / "resource.txt"
        resource.write_text("unrelated routing fixture\n")
        self.config = {key: str(resource) for key in adapter.PRIORITY_KEYS}
        self.config.update(assembly="hg19", vep_installed_version="116", clinvar_gene_stat="baseline.pkl.gz",
                           repeat_region_file_name="repeat.bed.gz", dispensable_gene_list=str(resource),
                           vep_cache_dir=str(self.folder / "cache"),
                           gene_mechanism_builder_script=str(self.source / "scripts/builder.py"),
                           interpro_mapping_pickle=str(resource), loeuf_prescore=str(resource))
        for relative in ("data/hpo/genes_to_phenotype.assertions.tsv.gz",
                         "data/alphamissense/alphamissense_mean_score.tsv",
                         "data/repeats/repeat.bed.gz", "data/hgnc/non_alt_loci_set.tsv"):
            path = self.source / relative
            path.parent.mkdir(parents=True, exist_ok=True)
            path.write_text("legacy fixture\n")
        self.baseline = self.source / "data/ClinVar/baseline.pkl.gz"
        self.selected = self.folder / "release_20261011/baseline.pkl.gz"
        self.write_stats(self.baseline, {"GENE": {"clinvar_pathogenic": 1}})
        self.counts = {"GENE": {"clinvar_pathogenic": 17}, "OTHER": {"clinvar_pathogenic": 3}}
        self.write_stats(self.selected, self.counts)

    @staticmethod
    def write_stats(path, counts):
        path.parent.mkdir(parents=True, exist_ok=True)
        with gzip.open(path, "wb") as handle:
            pickle.dump(counts, handle)

    def resources(self):
        return adapter.runtime_resources(self.config, self.source, "preannotated")

    def selected_record(self):
        paths = self.resources()
        manifest = adapter.inventory(paths, checksum=True)
        self.assertEqual(sum(key == "clinvar_gene_stat_used_by_shell" for key in manifest), 1)
        return manifest["clinvar_gene_stat_used_by_shell"]

    def test_manifest_records_actual_selected_file_identity_and_counts(self):
        before = self.baseline.read_bytes()
        self.config["clinvar_gene_stat"] = str(self.selected)
        record = self.selected_record()
        self.assertEqual(record["path"], str(self.selected))
        self.assertEqual(record["size"], self.selected.stat().st_size)
        self.assertEqual(record["sha256"], hashlib.sha256(self.selected.read_bytes()).hexdigest())
        with gzip.open(record["path"], "rb") as handle:
            self.assertEqual(pickle.load(handle), self.counts)
        self.assertEqual(self.baseline.read_bytes(), before)

    def test_caller_relative_directory_selection_has_same_absolute_manifest_identity(self):
        self.config["clinvar_gene_stat"] = "./" + str(self.selected.relative_to(self.folder))
        with working_directory(self.folder):
            record = self.selected_record()
        self.assertEqual(record["path"], str(self.selected))

    def test_bare_name_keeps_source_bundled_selection(self):
        self.assertEqual(self.selected_record()["path"], str(self.baseline))

    def test_explicit_missing_stat_is_not_replaced_by_same_named_baseline(self):
        self.config["clinvar_gene_stat"] = str(self.selected)
        self.selected.unlink()
        with self.assertRaisesRegex(adapter.ContractError, "clinvar_gene_stat_used_by_shell: missing"):
            self.resources()
        self.assertTrue(self.baseline.exists())

    def test_existing_nextflow_whitespace_rejection_is_preserved(self):
        selected = self.selected.with_name("selected gene statistics.pkl.gz")
        self.write_stats(selected, self.counts)
        self.config["clinvar_gene_stat"] = str(selected)
        with self.assertRaisesRegex(adapter.ContractError, "clinvar_gene_stat_used_by_shell:.*whitespace"):
            self.resources()

    def report(self, selected):
        self.config["clinvar_gene_stat"] = str(selected)
        config = self.folder / "profile.yaml"
        config.write_text(yaml.safe_dump(self.config))
        state = self.folder / "state"
        state.mkdir(exist_ok=True)
        output = self.folder / "readiness.json"
        command = ["check_cluster_readiness.py", "--config", str(config),
                   "--state-dir", str(state), "--output", str(output)]
        with patch.object(sys, "argv", command), contextlib.redirect_stdout(io.StringIO()):
            readiness.main()
        report = json.loads(output.read_text())
        records = [item for item in report["resources"] if item["key"] == "clinvar_gene_stat_used_by_shell"]
        self.assertEqual(len(records), 1)
        return records[0]

    def test_readiness_records_absolute_and_bare_selection(self):
        for value, expected in ((self.selected, self.selected), (self.baseline.name, self.baseline)):
            with self.subTest(value=value):
                record = self.report(value)
                self.assertEqual(record["path"], str(expected))
                self.assertEqual(record["status"], "present_not_scientifically_qualified")

    def test_readiness_resolves_caller_relative_directory_selection(self):
        with working_directory(self.folder):
            record = self.report(self.selected.relative_to(self.folder))
        self.assertEqual(record["path"], str(self.selected))
        self.assertEqual(record["status"], "present_not_scientifically_qualified")


if __name__ == "__main__":
    unittest.main()
