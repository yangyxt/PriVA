"""CPU allocations must not silently alter the legacy scientific profile."""
import contextlib
import importlib.util
import os
from pathlib import Path
import tempfile
import unittest

ROOT = Path(__file__).resolve().parents[2]
spec = importlib.util.spec_from_file_location("priva_cpu_adapter", ROOT / "bin/priva_task.py")
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


class CPUContract(unittest.TestCase):
    def setUp(self):
        temp = tempfile.TemporaryDirectory(dir=os.environ.get("PRIVA_TEST_TMPDIR"))
        self.addCleanup(temp.cleanup)
        self.folder = Path(temp.name)
        self.bundle = self.folder / "bundle"
        self.bundle.mkdir()
        (self.bundle / "family.ped").write_text("FAM\tchild\t0\t0\t0\t2\n")
        self.task = self.folder / "task"
        self.task.mkdir()

    def private(self, config, phase, cpus, mode="raw"):
        (self.bundle / "legacy.yaml").write_text(adapter.yaml.safe_dump(dict(assembly="hg38", **config)))
        original = (self.bundle / "legacy.yaml").read_bytes()
        contract = dict(base="input", assembly="hg38", mode=mode,
                        vcf_hub_mode="off", cadd_hub_mode="off")
        with working_directory(self.task):
            output, _ = adapter.private_config(self.bundle, contract, phase, cpus)
            result = adapter.read_config(output)
        self.assertEqual((self.bundle / "legacy.yaml").read_bytes(), original)
        return result

    def test_24_cpus_keep_qualified_22_forks_and_family_workers(self):
        result = self.private(dict(threads=22, threads_per_fam=4), "annotation", 24)
        self.assertEqual(result["threads"], 22)
        self.assertEqual(result["threads_per_fam"], 4)

    def test_explicit_string_worker_values_keep_their_types(self):
        result = self.private(dict(threads="22", threads_per_fam="4"), "annotation", 24)
        self.assertEqual((result["threads"], result["threads_per_fam"]), ("22", "4"))

    def test_insufficient_vep_allocation_fails_before_task_writes(self):
        with self.assertRaisesRegex(adapter.ContractError, "at least 24 allocated CPUs"):
            self.private(dict(threads=22, threads_per_fam=4), "annotation", 23)
        self.assertEqual(list(self.task.iterdir()), [])

    def test_larger_family_allocation_preserves_both_scientific_counts(self):
        result = self.private(dict(threads=22, threads_per_fam=4), "prioritization", 8)
        self.assertEqual((result["threads"], result["threads_per_fam"]), (22, 4))

    def test_insufficient_family_allocation_fails_before_task_writes(self):
        with self.assertRaisesRegex(adapter.ContractError, "at least 4 allocated CPUs"):
            self.private(dict(threads=22, threads_per_fam=4), "filtration", 3)
        self.assertEqual(list(self.task.iterdir()), [])

    def test_annotated_import_needs_no_vep_workers(self):
        result = self.private(dict(threads=22, threads_per_fam=4), "annotation", 1, mode="annotated")
        self.assertEqual((result["threads"], result["threads_per_fam"]), (22, 4))

    def test_raw_annotation_requires_explicit_scientific_worker_count(self):
        with self.assertRaisesRegex(adapter.ContractError, "explicit legacy-config threads"):
            self.private(dict(threads_per_fam=4), "annotation", 24)
        self.assertEqual(list(self.task.iterdir()), [])

    def test_invalid_worker_and_allocation_values_are_rejected(self):
        for value in (True, False, 0, -1, None, "", "2.0", "auto"):
            with self.subTest(worker=value), self.assertRaises(adapter.ContractError):
                adapter.task_thread_settings(dict(threads=value), "annotation", "raw", 24)
        for value in (True, 0, -1, 2.0, "24"):
            with self.subTest(cpus=value), self.assertRaises(adapter.ContractError):
                adapter.task_thread_settings(dict(threads=22), "annotation", "raw", value)


if __name__ == "__main__":
    unittest.main()
