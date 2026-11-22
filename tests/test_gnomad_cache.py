"""Real tiny indexed gnomAD families; no publisher downloads or score changes."""
import argparse
import importlib.util
import json
import os
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest

import pysam
import yaml

import fake_priva_binary


ROOT = Path(__file__).resolve().parents[1]
SPEC = importlib.util.spec_from_file_location("gnomad_family_tests", ROOT / "scripts/gnomad_cache.py")
CACHE = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(CACHE)
SPEC = importlib.util.spec_from_file_location("gnomad_adoption_tests", ROOT / "scripts/installer_config_paths.py")
PATHS = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(PATHS)


class GnomadFamilyTests(unittest.TestCase):
    def setUp(self):
        temporary = tempfile.TemporaryDirectory(prefix="priva_gnomad_family_", dir=os.environ.get("PRIVA_TEST_TMPDIR"))
        self.addCleanup(temporary.cleanup)
        self.folder = Path(temporary.name)
        self.source = self.folder / "selected_source"
        self.source.mkdir()
        self.managed = self.folder / "managed"
        self.reference = self.folder / "reference.fa"
        self.reference.write_text(">chr1\nA\n")
        self.config = self.folder / "config.json"
        self.settings = dict(cache_root=str(self.managed), tmp_dir=str(self.folder), assembly="hg38",
                             vep_installed_version="116", ref_genome=str(self.reference),
                             gnomad_vcf_dir=str(self.managed / "gnomAD"))
        self.write_family("hg38")

    def write_member(self, chromosome, *, missing_info=None, length=None):
        path = self.family[chromosome]
        header = "##fileformat=VCFv4.2\n" + "".join(
            f"##contig=<ID={name},length={length if name == chromosome and length is not None else size}>\n"
            for name, size in self.lengths.items())
        for name in CACHE.REQUIRED_INFO:
            if name != missing_info:
                kind = "Float" if name.startswith("AF_") else "Integer"
                header += f'##INFO=<ID={name},Number=1,Type={kind},Description="fixture">\n'
        plain = Path(str(path) + ".uncompressed.vcf")
        plain.write_text(header + "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n" +
                         f"{chromosome}\t1\t.\tA\tC\t.\tPASS\t.\n")
        pysam.tabix_compress(str(plain), str(path), force=True)
        pysam.tabix_index(str(path), preset="vcf", force=True)
        plain.unlink()

    def write_family(self, assembly):
        self.lengths = {chromosome: 100000 + number for number, chromosome in enumerate(CACHE.CHROMOSOMES)}
        self.lengths["chr1"] = CACHE.CHR1_LENGTH[assembly]
        Path(str(self.reference) + ".fai").write_text("".join(
            f"{name}\t{size}\t6\t1\t2\n" for name, size in self.lengths.items()))
        self.family = CACHE.family_paths(self.source / "selected.chrX.vcf.bgz")
        for chromosome in CACHE.CHROMOSOMES:
            self.write_member(chromosome)
        self.settings.update(assembly=assembly, gnomad_vcf_chrX=str(self.family["chrX"]))
        self.config.write_text(json.dumps(self.settings))

    def validate(self, assembly=None):
        return CACHE.validate_family(chr_x=self.family["chrX"], assembly=assembly or self.settings["assembly"], reference=self.reference)

    def snapshot(self):
        return {path: (CACHE.identity(path), path.read_bytes()) for path in self.source.iterdir()}


    def test_complete_family_checks_headers_indexes_and_leaves_all_source_bytes_unchanged(self):
        before = self.snapshot()
        result = self.validate()
        self.assertEqual(len(result["members"]), 24)
        self.assertFalse(result["scientific_acceptance"])
        self.assertEqual(self.snapshot(), before)


    def test_missing_required_joint_field_is_rejected(self):
        self.write_member("chr2", missing_info="AF_joint_nfe_XX")
        with self.assertRaisesRegex(ValueError, "caller-required INFO.*AF_joint_nfe_XX"):
            self.validate()

    def test_wrong_assembly_or_member_length_is_rejected(self):
        with self.assertRaisesRegex(ValueError, "requested gnomAD assembly|requested assembly"):
            self.validate("hg19")
        self.write_member("chr2", length=1)
        with self.assertRaisesRegex(ValueError, "differs from the selected assembly"):
            self.validate()

    def test_corrupt_index_is_rejected_without_repairing_originals(self):
        index = Path(str(self.family["chr2"]) + ".tbi")
        index.write_bytes(b"invalid index")
        before = self.snapshot()
        with self.assertRaises(ValueError):
            self.validate()
        self.assertEqual(self.snapshot(), before)

    def test_planner_adopts_only_the_exact_complete_family_and_known_companions(self):
        unrelated = self.source / "unrelated_payload.bin"
        unrelated.write_bytes(b"must not be copied")
        before = self.snapshot()
        PATHS.resolve_config(argparse.Namespace(config=self.config, source_root=ROOT, check=False, output=None, managed_updates_json=None))
        configured = yaml.safe_load(self.config.read_text())
        target = CACHE.family_paths(configured["gnomad_vcf_chrX"])
        self.assertTrue(target["chrX"].is_relative_to(self.managed))
        for chromosome, source in self.family.items():
            for suffix in ("", ".tbi"):
                original, copied = Path(str(source) + suffix), Path(str(target[chromosome]) + suffix)
                self.assertEqual(copied.read_bytes(), original.read_bytes())
                self.assertNotEqual(copied.stat().st_ino, original.stat().st_ino)
        self.assertFalse((target["chrX"].parent / unrelated.name).exists())
        self.assertEqual(self.snapshot(), before)

    def test_partial_custom_family_fails_before_planner_writes_or_creates_cache(self):
        self.family["chrY"].unlink()
        before = self.config.read_bytes()
        with self.assertRaisesRegex(ValueError, "member is missing"):
            PATHS.resolve_config(argparse.Namespace(config=self.config, source_root=ROOT, check=False, output=None, managed_updates_json=None))
        self.assertEqual(self.config.read_bytes(), before)
        self.assertFalse(self.managed.exists())


if __name__ == "__main__":
    unittest.main()
