"""Original-input and call contracts; fixture-only, without resource downloads."""
import gzip
import hashlib
import importlib.util
import json
import os
from pathlib import Path
import pickle
import shutil
import subprocess
import sys
import tempfile
import unittest
from unittest import mock

import fake_priva_binary


ROOT = Path(__file__).resolve().parents[1]
HELPER = ROOT / "scripts/installer_prerequisites.py"
SPEC = importlib.util.spec_from_file_location("installer_prerequisites", HELPER)
PREREQS = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(PREREQS)


class InstallerPrerequisites(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory(prefix="priva_prerequisites_",
                                               dir=os.environ.get("PRIVA_TEST_TMPDIR", os.environ.get("TMPDIR")))
        self.addCleanup(self.temp.cleanup)
        self.folder = Path(self.temp.name)
        self.source = self.folder / "selected_source"
        self.source.mkdir()
        self.mapping = self.folder / "original mapping.tsv"
        self.mapping.write_text(
            "Gene stable ID\tTranscript stable ID\tHGNC symbol\tUniProtKB/Swiss-Prot ID\tRefSeq mRNA ID\n"
            "ENSG000001\tENST000001.4\tGENE\tP00001\tNM_000001.2\n"
        )
        self.references = {}
        for assembly, length in PREREQS.CHR1_LENGTH.items():
            path = self.folder / f"original.{assembly}.fa"
            path.write_text(">chr1\nA\n")
            Path(str(path) + ".fai").write_text(f"chr1\t{length}\t6\t1\t2\n")
            self.references[assembly] = path
        self.chain = self.folder / "hg19ToHg38.over.chain"
        self.chain.write_text("chain 1 chr1 249250621 + 0 1 chr1 248956422 + 0 1 1\n1\n")
        self.cache = self.folder / "VEP_cache"
        for assembly in ("GRCh37", "GRCh38"):
            for relative in ("info.txt", "chr_synonyms.txt", "1/all_vars.gz", "1/all_vars.gz.csi"):
                path = self.cache / "homo_sapiens_merged" / f"113_{assembly}" / relative
                path.parent.mkdir(parents=True, exist_ok=True)
                path.write_text(f"species\thomo_sapiens\nassembly\t{assembly}\nvar_type\ttabix\n" if relative == "info.txt" else "fixture\n")
        self.stats = self.source / "data/ClinVar/clinvar_2star_stats.pkl.gz"
        self.stats.parent.mkdir(parents=True)
        with gzip.open(self.stats, "wb") as handle:
            pickle.dump({"GENE": {"ensembl_id": "ENSG1", "clinvar_pathogenic": set(),
                                 "large_aachange": set(), "small_aachange": set(),
                                 "putative_nmd_variants": set()}}, handle)
        self.hgnc = self.folder / "hgnc.tsv"
        self.hgnc.write_text("fixture\n")
        self.config = {
            "assembly": "hg38", "ref_genome": str(self.references["hg38"]),
            "nm_to_enst_tsv": str(self.mapping), "vep_installed_version": "113",
            "gofcards_hg19_fasta": str(self.references["hg19"]),
            "gofcards_hg38_fasta": str(self.references["hg38"]),
            "gofcards_hg19_to_hg38_chain": str(self.chain),
            "vep_cache_dir": str(self.cache), "gofcards_workdir": str(self.folder / "build_work"),
            "hgnc_table": str(self.hgnc), "threads": 2,
            "alphamissense_intolerant_motifs": str(self.folder / "motifs.pkl"),
            "alphamissense_vep_vcf": str(self.folder / "input.vcf.gz"),
            "alphamissense_pd_stat": str(self.folder / "input.stats.pkl"),
        }
        self.config_file = self.folder / "config.json"
        self.env = dict(os.environ, TEST_CONFIG=str(self.config_file), TEST_FOLDER=str(self.folder),
                        TEST_PYTHON=sys.executable, PYTHONDONTWRITEBYTECODE="1",
                        PRIVA_INSTALL_INSPECT_ROOTS=os.pathsep.join((str(ROOT), str(self.folder))),
                        # Builders run as main_install runs them inside its private
                        # view; native .rkyv values would otherwise start one.
                        PRIVA_CACHE_INSTALL_INTERNAL="1", **fake_priva_binary.environment(self.folder))
        self.env.pop("PRIVA_SMOKE_SKIP_VEP", None)


    def test_explicit_hpo_output_can_be_data_side_without_creating_or_replacing_it(self):
        target = self.folder / "data side's HPO" / "assertions.tsv.gz"
        self.config["hpo_assertions"] = str(target)
        PREREQS.validate_main(self.config, self.source)
        self.assertFalse(target.exists())  # Generated later by the original builder.
        target.parent.mkdir()
        target.write_bytes(b"explicit existing source remains unchanged\n")
        PREREQS.validate_main(self.config, self.source)
        self.assertEqual(target.read_bytes(), b"explicit existing source remains unchanged\n")

    def test_unset_hpo_destination_retains_legacy_default_without_creating_files(self):
        for value in (None, "", "null"):
            with self.subTest(value=value):
                self.config["hpo_assertions"] = value
                PREREQS.validate_main(self.config, self.source)
                self.assertFalse((self.source / "data/hpo").exists())

    def test_hpo_destination_keeps_old_host_type_and_inspection_guards(self):
        for value in ("/paedyl01/original/assertions.tsv.gz", True, False):
            with self.subTest(value=value):
                self.config["hpo_assertions"] = value
                with self.assertRaises(ValueError):
                    PREREQS.validate_main(self.config, self.source)
        self.config["hpo_assertions"] = str(self.folder / "data-side/assertions.tsv.gz")
        with mock.patch.dict(os.environ, PRIVA_INSTALL_INSPECT_ROOTS=os.pathsep.join((str(self.source), str(self.mapping)))), \
                mock.patch.object(PREREQS, "validate_mapping"), \
                self.assertRaisesRegex(ValueError, "outside permitted inspection roots.*hpo_assertions"):
            PREREQS.validate_main(self.config, self.source)

    def test_old_host_paths_are_reported_without_opening_them(self):
        self.config["gofcards_hg38_fasta"] = "/paedyl01/original/hg38.fa"
        with mock.patch.object(PREREQS, "validate_nm_to_enst", side_effect=AssertionError("opened data")):
            with self.assertRaisesRegex(ValueError, "Unresolved old-host path gofcards_hg38_fasta"):
                PREREQS.validate_main(self.config, self.source)

    def test_explicit_inspection_policy_rejects_paths_and_escaping_symlinks(self):
        with mock.patch.dict(os.environ, {"PRIVA_INSTALL_INSPECT_ROOTS": str(self.folder)}):
            with mock.patch.object(Path, "is_file", side_effect=AssertionError("inspected forbidden path")):
                with self.assertRaisesRegex(ValueError, "outside permitted"):
                    PREREQS.readable_file("/unauthorized/input.tsv", "nm_to_enst_tsv")
            link = self.folder / "outside"
            link.symlink_to("/unauthorized/input.tsv")
            with self.assertRaisesRegex(ValueError, "Symlink leaves"):
                PREREQS.readable_file(str(link), "nm_to_enst_tsv")

    def test_mapping_validation_preserves_original_bytes(self):
        before = self.mapping.read_bytes()
        PREREQS.validate_nm_to_enst(str(self.mapping))
        self.assertEqual(self.mapping.read_bytes(), before)
        self.mapping.write_text("RefSeq\tEnsembl\nNM_1\tENST1\n")
        with self.assertRaisesRegex(ValueError, "original five-column"):
            PREREQS.validate_nm_to_enst(str(self.mapping))

    def test_missing_bundled_clinvar_stats_reports_baseline_and_explicit_candidate_route(self):
        PREREQS.validate_main(self.config, self.source)
        self.stats.unlink()
        with self.assertRaisesRegex(ValueError, "ClinVar_Gene_stat_regenerate CONFIG SUMMARY MANE OUTPUT") as error:
            PREREQS.validate_main(self.config, self.source)
        self.assertIn("explicit ClinVar variant_summary and MANE inputs", str(error.exception))
        self.assertIn("does not replace the baseline", str(error.exception))
        self.assertNotIn("hardcoded old-host", str(error.exception))

    def test_snapshot_preflight_uses_only_existence_and_nonzero_size(self):
        # Intentionally opaque bytes: parsing/compressed integrity is outside
        # the user-requested policy for original repository-bundled snapshots.
        selected = {"assembly": "hg38"}
        files = [
            "HCSeeker/HC_spots.hg38.tsv", "gnomAD_RMC/gnomad_v2.1.1_rmc_flat.hg38.tsv",
            "alphamissense/alphamissense_mean_score.tsv", "vep/csq_interpretations.tsv",
            "repeats/ucsc.hg38.fasta.repeat.bed.gz", "hgnc/hgnc_non_alt_loci_set_sha256_890b18743cc8.tsv",
            "dispensable_genes/dispensable_gene_list.txt", "gofcards/gofcards_mechanism_reviews.tsv",
            "mondo/disease_scope_overrides.tsv",
            "patho_mechanism/gene_nonlof_mechanism_curated_assertions.schema.json",
        ]
        for relative in files:
            path = self.source / "data" / relative
            path.parent.mkdir(parents=True, exist_ok=True)
            path.write_bytes(b"opaque original snapshot")
        PREREQS.validate_bundled_snapshots(selected, self.source)
        last = self.source / "data" / files[-1]
        last.write_bytes(b"")
        with self.assertRaisesRegex(ValueError, "Missing, empty"):
            PREREQS.validate_bundled_snapshots(selected, self.source)
        last.unlink()
        with self.assertRaisesRegex(ValueError, "Missing, empty"):
            PREREQS.validate_bundled_snapshots(selected, self.source)

    def test_install_preflight_defers_only_absent_or_empty_genuine_bundled_baseline(self):
        for state in ("missing", "empty"):
            with self.subTest(state=state):
                if self.stats.exists():
                    self.stats.unlink()
                if state == "empty":
                    self.stats.touch()
                PREREQS.validate_main(self.config, self.source, allow_missing_bundled_clinvar=True)
                with self.assertRaises(ValueError):
                    PREREQS.validate_main(self.config, self.source)
        external = self.folder / "custom-missing.pkl.gz"
        with self.assertRaises((ValueError, OSError)):
            PREREQS.validate_main(dict(self.config, clinvar_gene_stat=str(external)),
                                  self.source, allow_missing_bundled_clinvar=True)
        self.stats.unlink()
        self.stats.symlink_to(external)
        with self.assertRaises((ValueError, OSError)):
            PREREQS.validate_main(self.config, self.source, allow_missing_bundled_clinvar=True)
        self.assertFalse(external.exists())


    def test_explicit_generated_clinvar_stat_is_checked_at_its_selected_absolute_path(self):
        selected = self.folder / "generated statistics" / "release.pkl.gz"
        selected.parent.mkdir()
        selected.write_bytes(self.stats.read_bytes())
        self.config["clinvar_gene_stat"] = str(selected)
        self.assertEqual(PREREQS.clinvar_stats_path(self.config, self.source), selected)
        PREREQS.validate_selected_clinvar_stats(self.config, self.source)
        selected.unlink()
        with self.assertRaisesRegex(ValueError, "Cannot reuse bundled ClinVar gene statistics"):
            PREREQS.validate_selected_clinvar_stats(self.config, self.source)

    def test_generated_native_clinvar_stat_archive_is_checked_by_its_header(self):
        # priva clinvar-gene-stats writes a PRIVARKV archive with schema 5.
        selected = self.folder / "generated statistics" / "release.rkyv"
        selected.parent.mkdir()
        archive = b"\0" * 16
        def header(schema, length=len(archive)):
            return (b"PRIVARKV" + (1).to_bytes(4, "little") + schema.to_bytes(4, "little")
                    + (1).to_bytes(4, "little") + (1).to_bytes(4, "little") + length.to_bytes(8, "little"))
        selected.write_bytes(header(5) + archive)
        self.config["clinvar_gene_stat"] = str(selected)
        PREREQS.validate_selected_clinvar_stats(self.config, self.source)
        for value in (header(4) + archive, header(5, len(archive) + 1) + archive, self.stats.read_bytes()):
            selected.write_bytes(value)
            with self.assertRaisesRegex(ValueError, "Cannot reuse bundled ClinVar gene statistics"):
                PREREQS.validate_selected_clinvar_stats(self.config, self.source)

    def test_clinvar_stat_relative_directory_is_retained_but_bare_name_is_bundled_relative(self):
        self.config["clinvar_gene_stat"] = "generated/release.pkl.gz"
        self.assertEqual(PREREQS.clinvar_stats_path(self.config, self.source), Path("generated/release.pkl.gz"))
        self.config["clinvar_gene_stat"] = "release.pkl.gz"
        self.assertEqual(PREREQS.clinvar_stats_path(self.config, self.source), self.source / "data/ClinVar/release.pkl.gz")

    def test_both_merged_caches_and_correct_assemblies_are_required_for_build(self):
        PREREQS.validate_gofcards(self.config)
        self.config["gofcards_hg19_fasta"] = self.config["ref_genome"]
        with self.assertRaisesRegex(ValueError, "hg19 chromosome 1"):
            PREREQS.validate_gofcards(self.config)
        self.config["gofcards_hg19_fasta"] = str(self.references["hg19"])
        (self.cache / "homo_sapiens_merged/113_GRCh37/1/all_vars.gz").unlink()
        with self.assertRaisesRegex(ValueError, "GRCh37 merged VEP 113"):
            PREREQS.validate_gofcards(self.config)

    def test_reverse_chain_is_rejected(self):
        self.chain.write_text("chain 1 chr1 248956422 + 0 1 chr1 249250621 + 0 1 1\n1\n")
        with self.assertRaisesRegex(ValueError, "forward UCSC"):
            PREREQS.validate_gofcards(self.config)

    def test_main_allows_selected_reference_and_cache_to_be_installed(self):
        selected = self.references["hg38"]
        selected.unlink()
        Path(str(selected) + ".fai").unlink()
        shutil.rmtree(self.cache / "homo_sapiens_merged/113_GRCh38")
        PREREQS.validate_gofcards(self.config, main_preflight=True, allow_missing_merged_caches=True)
        with self.assertRaisesRegex(ValueError, "gofcards_hg38_fasta"):
            PREREQS.validate_gofcards(self.config)

    def test_main_rejects_existing_wrong_assembly_before_costly_work(self):
        self.config["ref_genome"] = str(self.references["hg19"])
        self.config["gofcards_hg38_fasta"] = self.config["ref_genome"]
        with self.assertRaisesRegex(ValueError, "hg38 chromosome 1"):
            PREREQS.validate_gofcards(self.config, main_preflight=True)


    def test_compatible_bundle_default_is_resolved_without_config_mutation(self):
        self.config["assembly"] = "hg19"
        self.config.pop("nm_to_enst_tsv")
        before = dict(self.config)
        resolved = PREREQS.validate_mapping(self.config, ROOT)
        self.assertEqual(resolved, ROOT / PREREQS.BIOMART.RELATIVE_PATH)
        self.assertEqual(self.config, before)


    def test_missing_bundle_can_be_deferred_only_during_main_preflight(self):
        self.config.pop("nm_to_enst_tsv")
        for assembly in ("hg19", "hg38"):
            with self.subTest(assembly=assembly):
                self.config["assembly"] = assembly
                PREREQS.validate_main(self.config, self.source, allow_missing_bundled_mapping=True)
                with self.assertRaisesRegex(ValueError, "nm_to_enst_tsv"):
                    PREREQS.validate_main(self.config, self.source)

    def test_grch38_bundle_default_is_verified_offline_without_config_mutation(self):
        self.config.pop("nm_to_enst_tsv")
        for assembly in ("hg38", "GRCh38"):
            with self.subTest(assembly=assembly):
                self.config["assembly"] = assembly
                before = dict(self.config)
                source = PREREQS.BIOMART.source_for_assembly(assembly)
                with mock.patch.object(PREREQS.BIOMART, "urlopen", side_effect=AssertionError("unexpected network")):
                    resolved = PREREQS.validate_mapping(self.config, ROOT)
                self.assertEqual(resolved, ROOT / source["relative_path"])
                self.assertEqual(source["source_id"], "ensembl-grch38-mart113")
                self.assertEqual(self.config, before)


    def test_grch38_pinned_basename_requires_its_provenance(self):
        source = PREREQS.BIOMART.source_for_assembly("hg38")
        path = self.folder / source["basename"]
        path.write_bytes(self.mapping.read_bytes())
        with self.assertRaises(OSError):
            PREREQS.validate_nm_to_enst(str(path), "hg38")

    def test_missing_default_symlink_never_becomes_an_original_download_target(self):
        self.config.pop("nm_to_enst_tsv")
        source = PREREQS.BIOMART.source_for_assembly("hg38")
        bundled = self.source / source["relative_path"]
        bundled.parent.mkdir(parents=True, exist_ok=True)
        original = self.folder / "missing_original.tsv"
        bundled.symlink_to(original)
        with mock.patch.dict(os.environ, {"PRIVA_INSTALL_INSPECT_ROOTS": str(self.folder)}):
            with self.assertRaisesRegex(ValueError, "Explicit nm_to_enst_tsv is missing"):
                PREREQS.mapping_fetch_source(self.config, self.source)
        self.assertFalse(original.exists())
        self.assertTrue(bundled.is_symlink())


if __name__ == "__main__":
    unittest.main()
