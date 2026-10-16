import contextlib
import gzip
import hashlib
import importlib.util
import json
import os
from pathlib import Path
import tempfile
import unittest
from types import SimpleNamespace

ROOT = Path(__file__).resolve().parents[2]
spec = importlib.util.spec_from_file_location("priva_task", ROOT / "bin/priva_task.py")
adapter = importlib.util.module_from_spec(spec)
spec.loader.exec_module(adapter)


@contextlib.contextmanager
def working_directory(path):
    original = Path.cwd()
    os.chdir(path)
    try:
        yield
    finally:
        os.chdir(original)


class AdapterContracts(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory(dir=ROOT / "tests/nextflow")
        self.addCleanup(self.temp.cleanup)
        self.folder = Path(self.temp.name)

    def test_source_digest_uses_nextflow_path_string_order(self):
        source = self.folder / "source-order"
        (source / "conf/cadd").mkdir(parents=True)
        for name, content in (("conf/cadd.config", "flat\n"),
                              ("conf/cadd/profile.json", "nested\n"),
                              ("main.nf", "main\n"),
                              ("nextflow.config", "settings\n")):
            (source / name).write_text(content)
        # Explicit lexical path order, including the dot/slash prefix collision
        # that made the Java and Python provenance digests disagree in CI.
        expected = hashlib.sha256(
            b"conf/cadd.config\0flat\nconf/cadd/profile.json\0nested\n"
            b"main.nf\0main\nnextflow.config\0settings\n").hexdigest()
        self.assertEqual(adapter.source_digest(source), expected)

    def test_ped_preserves_cohort_context_and_selects_family(self):
        path = self.folder / "cohort.ped"
        text = "OTHER\tother\t0\t0\t2\t2\nFAM\tchild\tdad\tmom\t1\t2\tHP:0000001\nFAM\tmom\t0\t0\t2\t1\nFAM\tdad\t0\t0\t1\t1\n"
        path.write_text(text)
        self.assertEqual(adapter.validate_ped(path, "FAM", ["child", "dad", "mom", "other"]), ["child", "mom", "dad"])
        self.assertEqual(path.read_text(), text)

    def test_ped_rejects_missing_family_samples_unknown_phenotype_and_wrong_proband(self):
        path = self.folder / "family.ped"
        for contents, samples in [("FAM\tchild\t0\t0\t1\t2\n", []),
                                  ("OTHER\tchild\t0\t0\t1\t2\n", ["child"]),
                                  ("FAM\tchild\t0\t0\t1\t0\n", ["child"]),
                                  ("FAM\tchild\t0\t0\t1\t1\n", ["child"])]:
            path.write_text(contents)
            with self.assertRaises(adapter.ContractError):
                adapter.validate_ped(path, "FAM", samples)

    def test_duplicate_yaml_is_not_silently_overwritten(self):
        path = self.folder / "bad.yaml"
        path.write_text("assembly: hg19\nassembly: hg38\n")
        with self.assertRaises(adapter.ContractError):
            adapter.read_config(path)

    def test_missing_and_relative_cache_paths_fail(self):
        for value in ("missing/cache", str(self.folder / "missing"), "", "/path with space/cache"):
            with self.assertRaises(adapter.ContractError):
                adapter.checked_path(value, "cache")

    def test_lfs_pointer_is_rejected(self):
        path = self.folder / "pointer.tsv"
        path.write_text("version https://git-lfs.github.com/spec/v1\noid sha256:abc\nsize 123\n")
        with self.assertRaises(adapter.ContractError):
            adapter.checked_path(str(path), "cache")

    def test_copy_is_private_even_when_original_is_read_only(self):
        source = self.folder / "source.tsv"
        source.write_text("original\n")
        source.chmod(0o444)
        copied = adapter.copy_writable(source, self.folder / "task/copy.tsv")
        copied.write_text("changed\n")
        self.assertEqual(source.read_text(), "original\n")
        self.assertNotEqual(source.stat().st_ino, copied.stat().st_ino)

    def test_inventory_tracks_nested_cache_changes(self):
        root = self.folder / "cache"
        root.mkdir()
        (root / "release").mkdir()
        path = root / "release/data.tsv"
        path.write_text("before")
        before = adapter.inventory({"vep": root})
        path.write_text("after and different size")
        after = adapter.inventory({"vep": root})
        self.assertNotEqual(before["vep/release/data.tsv"], after["vep/release/data.tsv"])

    def test_cadd_header_contract(self):
        path = self.folder / "scores.tsv"
        path.write_text("#Chrom\tPos\tRef\tAlt\tRawScore\tPHRED\tFeatureID\n")
        adapter.validate_cadd(path)
        path.write_text("Chr\tStart\tRef\tAlt\n")
        with self.assertRaises(adapter.ContractError):
            adapter.validate_cadd(path)

    def test_private_config_keeps_scientific_values_and_isolates_mutable_resources(self):
        bundle = self.folder / "bundle"
        bundle.mkdir()
        (bundle / "family.ped").write_text("FAM\tchild\t0\t0\t1\t2\n")
        controls = self.folder / "controls.vcf.gz"
        controls.write_bytes(b"immutable-control-bytes")
        Path(str(controls) + ".tbi").write_bytes(b"index")
        cds = self.folder / "coding.fa.gz"
        cds.write_bytes(b"fasta")
        cache = self.folder / "coding.cds_cache.pkl.gz"
        cache.write_bytes(b"pickle")
        (bundle / "legacy.yaml").write_text(adapter.yaml.safe_dump({
            "assembly": "hg19", "af_cutoff": 0.05, "extreme_rare_PAF": 0.0001,
            "apply_clingen_override": False, "control_vcf": str(controls), "cds_fasta_file": str(cds),
            "clinvar_aa_stat": "/shared/same/cache.pkl", "hub_vcf_file": "/shared/hub.vcf.gz"}))
        contract = {"base": "input", "assembly": "hg19", "mode": "raw", "vcf_hub_mode": "off", "cadd_hub_mode": "off"}
        task = self.folder / "task"
        task.mkdir()
        with working_directory(task):
            config_path, temp = adapter.private_config(bundle, contract, "prioritization", 2)
            config = adapter.read_config(config_path)
        self.assertEqual(config["af_cutoff"], 0.05)
        self.assertFalse(config["apply_clingen_override"])
        self.assertEqual(config["clinvar_aa_stat"], "/shared/same/cache.pkl")
        self.assertEqual(config["threads_per_fam"], 2)
        self.assertEqual(config["hub_vcf_file"], "")
        self.assertTrue(str(config["control_vcf"]).startswith(str(task)))
        self.assertTrue(Path(config["cds_fasta_file"].replace(".fa.gz", ".cds_cache.pkl.gz")).exists())
        Path(config["control_vcf"]).write_bytes(b"modified")
        self.assertEqual(controls.read_bytes(), b"immutable-control-bytes")

    def test_legacy_success_without_outputs_fails_closed(self):
        with self.assertRaises(adapter.ContractError):
            adapter.require_outputs([self.folder / "missing.filtered.acmg.tsv"])

    def make_launch_contract(self):
        source = self.folder / "source"
        source.mkdir()
        for name in ("main.nf", "nextflow.config"):
            (source / name).write_text("frozen source\n")
        for relative in ("data/hpo/genes_to_phenotype.assertions.tsv.gz",
                         "data/alphamissense/alphamissense_mean_score.tsv",
                         "data/ClinVar/gene.pkl", "data/repeats/repeats.bed.gz",
                         "data/hgnc/non_alt_loci_set.tsv", "data/dispensable_genes/dispensable_gene_list.txt"):
            path = source / relative
            path.parent.mkdir(parents=True, exist_ok=True)
            path.write_bytes(b"synthetic resource fixture")
        cache = self.folder / "resource.pkl"
        cache.write_bytes(b"same fixed resource")
        config = {key: str(cache) for key in adapter.PRIORITY_KEYS}
        config.update(assembly="hg19", clinvar_gene_stat="gene.pkl", repeat_region_file_name="repeats.bed.gz")
        config_path = self.folder / "legacy.yaml"
        config_path.write_text(adapter.yaml.safe_dump(config))
        manifest_path = self.folder / "resources.json"
        adapter.build_manifest(SimpleNamespace(config=str(config_path), legacy_root=str(source), mode="annotated",
            vcf_hub_mode="off", cadd_hub_mode="copy", checksum=False, output=str(manifest_path)))
        ped = self.folder / "family.ped"
        ped.write_text("FAM\tchild\t0\t0\t1\t2\n")
        vcf = self.folder / "input.vcf.gz"
        vcf.write_bytes(b"fixture only: launch inventories inputs; validation parses VCF")
        index = self.folder / "input.vcf.gz.tbi"
        index.write_bytes(b"index")
        args = SimpleNamespace(config=str(config_path), legacy_root=str(source), mode="annotated",
            vcf_hub_mode="off", cadd_hub_mode="copy", resources=str(manifest_path),
            source_digest=adapter.source_digest(source), software_id="frozen-test-environment", family="FAM",
            assembly="hg19", vcf=str(vcf), vcf_index=str(index), ped=str(ped), cadd=None,
            annotation_cpus=1, family_cpus=1, stub=False, outdir=str(self.folder / "published"))
        return args, cache

    def test_launch_rechecks_resources_even_for_identical_resume(self):
        args, resource = self.make_launch_contract()
        adapter.launch(args)
        adapter.launch(args)  # Identical resume is allowed.
        resource.write_bytes(b"changed after completed cached run")
        with self.assertRaisesRegex(adapter.ContractError, "Resource changed"):
            adapter.launch(args)

    def test_output_directory_rejects_changed_family_identity(self):
        args, _ = self.make_launch_contract()
        adapter.launch(args)
        args.family = "SECOND_FAMILY"
        with self.assertRaisesRegex(adapter.ContractError, "different inputs"):
            adapter.launch(args)

    def test_source_changed_after_launch_is_rejected(self):
        args, _ = self.make_launch_contract()
        (Path(args.legacy_root) / "main.nf").write_text("changed source")
        with self.assertRaisesRegex(adapter.ContractError, "Source changed"):
            adapter.launch(args)

    def test_raw_scoring_gate_prevents_shared_conda_cache_mutation(self):
        with self.assertRaisesRegex(adapter.ContractError, "read-only"):
            adapter.require_qualified_execution("raw", False)
        adapter.require_qualified_execution("raw", True)
        adapter.require_qualified_execution("annotated", False)

    def test_space_delimited_ped_is_rejected_like_legacy_parser(self):
        path = self.folder / "family.ped"
        path.write_text("FAM child 0 0 1 2\n")
        with self.assertRaisesRegex(adapter.ContractError, "tab-delimited"):
            adapter.validate_ped(path, "FAM", ["child"])


if __name__ == "__main__":
    unittest.main()
