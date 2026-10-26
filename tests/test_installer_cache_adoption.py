"""Cache-root adoption preserves selected bytes; fixtures never use real caches."""
import argparse
import importlib.util
import json
import os
from pathlib import Path
import sys
import tempfile
import unittest
from unittest import mock

import yaml


ROOT = Path(__file__).resolve().parents[1]
SPEC = importlib.util.spec_from_file_location("installer_cache_adoption", ROOT / "scripts/installer_config_paths.py")
PATHS = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(PATHS)


class CacheAdoptionTests(unittest.TestCase):
    def setUp(self):
        temporary = tempfile.TemporaryDirectory(prefix="priva_cache_adoption_", dir=os.environ.get("PRIVA_TEST_TMPDIR"))
        self.addCleanup(temporary.cleanup)
        self.folder = Path(temporary.name)
        self.cache = self.folder / "new_cache"
        self.config = self.folder / "config.yaml"
        self.settings = dict(cache_root=str(self.cache), tmp_dir=str(self.folder / "temporary"),
                             assembly="hg38", vep_installed_version="116")

    @staticmethod
    def vcf_fixture(version=1, body=b""):
        return (f"##fileformat=VCFv4.2\n##reference=GRCh38\n##SpliceAIVersion={version}\n"
                "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n").encode() + body

    def write(self):
        self.config.write_text(yaml.safe_dump(self.settings))

    def resolve(self, *, check=False, output=None, source_root=ROOT):
        args = argparse.Namespace(config=self.config, source_root=source_root,
                                  managed_updates_json=None, check=check, output=output)
        PATHS.resolve_config(args)

    def test_old_plugin_root_is_adopted_with_index_and_provenance_before_config_rewrite(self):
        plugins = self.folder / "old_plugins"
        source = plugins / "AlphaMissense/AlphaMissense_hg38.tsv.gz"
        source.parent.mkdir(parents=True)
        source.write_bytes(b"selected original bytes; scientific validation is separate")
        Path(str(source) + ".tbi").write_bytes(b"selected original index")
        Path(str(source) + ".provenance.json").write_bytes(b'{"source":"unchanged"}\n')
        self.settings.update(vep_plugins_cachedir=str(plugins), alphamissense_prescore=str(source))
        self.write()
        before = PATHS.cache_identity(source)
        with mock.patch("builtins.print"):
            PATHS.configure_cache_paths(self.config, self.settings, "--check", ROOT)
        self.assertFalse(self.cache.exists())
        PATHS.configure_cache_paths(self.config, self.settings, "write", ROOT)
        configured = yaml.safe_load(self.config.read_text())
        target = Path(configured["alphamissense_prescore"])
        self.assertEqual(target, self.cache / "VEP_plugin_cache/AlphaMissense/supplied/alphamissense_prescore/AlphaMissense_hg38.tsv.gz")
        self.assertEqual(target.read_bytes(), source.read_bytes())
        self.assertNotEqual(target.stat().st_ino, source.stat().st_ino)
        self.assertEqual(PATHS.cache_identity(source), before)
        for suffix in (".tbi", ".provenance.json"):
            self.assertEqual(Path(str(target) + suffix).read_bytes(), Path(str(source) + suffix).read_bytes())
        receipt = json.loads((self.cache / "installer-source-adoptions.json").read_text())
        self.assertIn(str(source), [item["source"] for item in receipt["copies"]])

    def test_relocated_invalid_alpha_source_remains_explicit_and_cannot_auto_download(self):
        old = self.folder / "old_plugins"
        source = old / "AlphaMissense/AlphaMissense_hg38.tsv.gz"
        source.parent.mkdir(parents=True)
        source.write_bytes(b"invalid explicitly supplied source; never substitute it")
        self.settings.update(vep_plugins_cachedir=str(old), alphamissense_prescore=str(source))
        self.write()
        self.resolve()
        selected = yaml.safe_load(self.config.read_text())
        candidate = Path(selected["alphamissense_prescore"])
        canonical = Path(selected["vep_plugins_cachedir"]) / "AlphaMissense/AlphaMissense_hg38.tsv.gz"
        self.assertNotEqual(candidate, canonical)
        spec = importlib.util.spec_from_file_location("alpha_adoption_policy", ROOT / "scripts/alphamissense_cache.py")
        alpha = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(alpha)
        with mock.patch.object(alpha, "download", side_effect=AssertionError("must not download for an invalid explicit source")), \
                mock.patch.object(alpha.shutil, "which", return_value="fixture tool"):
            with self.assertRaisesRegex(ValueError, "publisher size mismatch"):
                alpha.install(cache_dir=Path(selected["vep_plugins_cachedir"]), assembly="hg38", source=candidate)
        self.assertEqual(candidate.read_bytes(), source.read_bytes())
        self.assertFalse((canonical.parent / "AlphaMissense_hg38.ready.json").exists())

    def test_explicit_custom_input_is_not_substituted_when_missing(self):
        missing = self.folder / "selected-release-custom.tsv.gz"
        self.settings["primateai_prescore"] = str(missing)
        self.write()
        before = self.config.read_bytes()
        with self.assertRaisesRegex(ValueError, "refusing to substitute another dataset/release"):
            self.resolve()
        self.assertEqual(self.config.read_bytes(), before)
        self.assertFalse(self.cache.exists())

    def test_custom_replacement_tables_preserve_selected_hgnc_and_missing_selections_fail(self):
        for key in sorted(PATHS.REPLACEMENT_TABLE_KEYS):
            with self.subTest(key=key):
                source = self.folder / (key + ".selected")
                source.write_bytes(b"selected replacement table; retain these exact bytes")
                self.settings[key] = str(source)
                self.write()
                self.resolve()
                target = Path(yaml.safe_load(self.config.read_text())[key])
                if key == "hgnc_table":
                    self.assertEqual(target, source)
                    self.assertEqual(target.stat().st_ino, source.stat().st_ino)
                    self.assertFalse((self.cache / "data/hgnc").exists())
                else:
                    self.assertTrue(target.is_relative_to(self.cache))
                    self.assertNotEqual(target.stat().st_ino, source.stat().st_ino)
                self.assertEqual(target.read_bytes(), source.read_bytes())
                source.unlink()
                self.write()
                before = self.config.read_bytes()
                # The GoFCards review is needed only when normalization must
                # rebuild. Force that route for this missing-source check.
                forced = {"PRIVA_FORCE_GOFCARDS_CACHE": "1"} if key == "gofcards_mechanism_review_tsv" else {}
                with mock.patch.dict(os.environ, forced):
                    with self.assertRaisesRegex(ValueError, "refusing to substitute another dataset/release"):
                        self.resolve()
                self.assertEqual(self.config.read_bytes(), before)
                del self.settings[key]

    def test_canonical_replacement_tables_need_only_existence_and_nonzero_size(self):
        source_root = self.folder / "table_checkout"
        for key in PATHS.REPLACEMENT_TABLE_KEYS:
            path = source_root / PATHS.BUNDLED_PATHS[key]
            path.parent.mkdir(parents=True, exist_ok=True)
            path.write_bytes(b"nonempty canonical bundled snapshot")
            self.settings[key] = PATHS.BUNDLED_PATHS[key]
        self.write()
        with mock.patch.object(PATHS, "sha256_file", side_effect=AssertionError("never hash bundles")), \
                mock.patch.object(PATHS.shutil, "copy2", side_effect=AssertionError("never copy bundles")):
            self.resolve(source_root=source_root)
        configured = yaml.safe_load(self.config.read_text())
        for key in PATHS.REPLACEMENT_TABLE_KEYS:
            self.assertEqual(configured[key], str(source_root / PATHS.BUNDLED_PATHS[key]))
        self.assertFalse(self.cache.exists())

    def test_existing_custom_basename_is_not_substituted_by_the_legacy_bundle_alias(self):
        source = self.folder / "Interpro_entry_mapping.pkl.gz"
        source.write_bytes(b"explicit selected replacement with the bundle basename")
        self.settings["interpro_mapping_pickle"] = source.name
        self.write()
        previous = Path.cwd()
        try:
            os.chdir(self.folder)
            self.resolve()
        finally:
            os.chdir(previous)
        target = Path(yaml.safe_load(self.config.read_text())["interpro_mapping_pickle"])
        self.assertTrue(target.is_relative_to(self.cache))
        self.assertEqual(target.read_bytes(), source.read_bytes())

    def test_missing_custom_table_basename_is_not_replaced_by_a_canonical_bundle(self):
        self.settings["hgnc_table"] = "non_alt_loci_set.tsv"
        self.write()
        previous = Path.cwd()
        try:
            os.chdir(self.folder)
            with self.assertRaisesRegex(ValueError, "refusing to substitute another dataset/release"):
                self.resolve()
        finally:
            os.chdir(previous)
        self.assertFalse(self.cache.exists())

    def test_missing_custom_runtime_inputs_already_under_cache_root_fail_before_writes(self):
        keys = ("primateai_prescore", "spliceai_snv_prescore", "spliceai_indel_prescore",
                "splicevault_prescore", "utr_annotator_file", "human_ancestor_fasta",
                "loftee_conservation_file", "gerp_bigwig", "loftee_repo", "cadd_base_dir",
                "mavedb_file", "cds_fasta_file", "conservation_file", "ref_genome",
                "gnomad_vcf_chrX", "nm_to_enst_tsv", "mondo_obo", "clinvar_vcv_xml")
        for key in keys:
            with self.subTest(key=key):
                self.settings[key] = str(self.cache / "selected_custom" / key)
                self.write()
                before = self.config.read_bytes()
                with self.assertRaisesRegex(ValueError, "refusing to substitute another dataset/release"):
                    self.resolve()
                self.assertEqual(self.config.read_bytes(), before)
                self.assertFalse(self.cache.exists())
                del self.settings[key]

    def test_exact_documented_missing_public_defaults_remain_acquirable(self):
        defaults = PATHS.portable_updates(self.settings, ROOT)
        for key in PATHS.AUTO_SOURCE_KEYS:
            if key in defaults:
                self.settings[key] = defaults[key]
        self.write()
        before = self.config.read_bytes()
        with mock.patch("builtins.print"):
            self.resolve(check=True)
        self.assertEqual(self.config.read_bytes(), before)
        self.assertFalse(self.cache.exists())

    def test_disabled_pext_missing_sources_are_dormant(self):
        self.settings.update(pext_enabled=False, pext_mean_bw=str(self.folder / "missing_selected.bw"),
                             pext_tissue_bw_dir=str(self.folder / "missing_selected_tissues"))
        self.write()
        self.resolve()
        configured = yaml.safe_load(self.config.read_text())
        self.assertEqual(configured["pext_mean_bw"], self.settings["pext_mean_bw"])
        self.assertEqual(configured["pext_tissue_bw_dir"], self.settings["pext_tissue_bw_dir"])
        self.assertFalse(self.cache.exists())

    def test_missing_custom_gofcards_review_is_deferred_and_existing_workdir_stays_visible_in_preflight(self):
        work = self.folder / "original_build_work"
        work.mkdir()
        (work / "gofcards_exact_gof.normalized.json.gz").write_bytes(b"opaque selected normalized cache")
        review = self.folder / "unused_missing_review.tsv"
        self.settings.update(gofcards_workdir=str(work), gofcards_mechanism_review_tsv=str(review))
        self.write()
        planning = self.folder / "planning.yaml"
        self.resolve(output=planning)
        prepared = yaml.safe_load(planning.read_text())
        self.assertEqual(prepared["gofcards_workdir"], str(work))
        self.assertEqual(prepared["gofcards_mechanism_review_tsv"], str(review))
        self.assertFalse(self.cache.exists())
        self.resolve()
        configured = yaml.safe_load(self.config.read_text())
        self.assertTrue(Path(configured["gofcards_workdir"]).is_relative_to(self.cache))
        self.assertEqual(configured["gofcards_mechanism_review_tsv"], str(review))
        self.assertFalse(review.exists())

    def test_missing_or_empty_default_clinvar_survives_planning_preflight_and_refresh_decision(self):
        spec = importlib.util.spec_from_file_location("clinvar_adoption_policy", ROOT / "scripts/clinvar_gene_refresh.py")
        refresh = importlib.util.module_from_spec(spec)
        with mock.patch.object(sys, "path", [str(ROOT / "scripts"), *sys.path]):
            spec.loader.exec_module(refresh)
        source_root = self.folder / "fresh_checkout"
        baseline = source_root / "data/ClinVar/clinvar_2star_stats.pkl.gz"
        baseline.parent.mkdir(parents=True)
        registry = {"schema_version": 1, "resources": {
            "clinvar_gene_statistics": {"refresh": {"threshold_months": 3, "automatic": True,
                "policy": "latest_filedate_vs_bundled_mtime"}},
            "clinvar_vcf": {"source": {"latest_url_template": refresh.NCBI + "/pub/clinvar/vcf_{assembly}/clinvar.vcf.gz"}},
            "clinvar_variant_summary": {"source": {"url": refresh.SUMMARY_URL}},
            "mane_select": {"source": {"url": refresh.MANE_URL}},
        }}
        (source_root / "cache_resources.json").write_text(json.dumps(registry))
        mapping = self.folder / "selected_mapping.tsv"
        mapping.write_text("Gene stable ID\tTranscript stable ID\tProtein stable ID\tChromosome/scaffold name\tRefSeq mRNA ID\n"
                           "ENSG1\tENST1\tENSP1\t1\tNM_1\n")
        self.settings["nm_to_enst_tsv"] = str(mapping)
        latest = {"file_date": "2026-04-30", "publisher_md5": "1" * 32}
        for state in ("missing", "empty"):
            for selected in ("", baseline.name, str(baseline)):
                with self.subTest(state=state, selected=selected):
                    baseline.unlink(missing_ok=True)
                    if state == "empty":
                        baseline.touch()
                    self.settings["clinvar_gene_stat"] = selected
                    self.write()
                    original = self.config.read_bytes()
                    planning = self.folder / "planning.yaml"
                    self.resolve(output=planning, source_root=source_root)
                    self.assertEqual(self.config.read_bytes(), original)
                    self.assertEqual(yaml.safe_load(planning.read_text())["clinvar_gene_stat"], str(baseline))
                    for config in (planning, self.config):
                        if config == self.config:
                            self.resolve(source_root=source_root)
                        prepared = yaml.safe_load(config.read_text())
                        self.assertEqual(prepared["clinvar_gene_stat"], str(baseline))
                        refresh.validate_main(prepared, source_root, allow_missing_bundled_clinvar=True)
                        with mock.patch.object(refresh, "latest_metadata", return_value=latest):
                            plan = refresh.plan_refresh(config, source_root)
                        self.assertEqual(plan["action"], "rebuild")
                        self.assertTrue(plan["baseline_missing"])
                        self.assertTrue(Path(plan["output_path"]).is_relative_to(self.cache))
                    self.assertEqual(baseline.exists(), state == "empty")
                    if baseline.exists():
                        self.assertEqual(baseline.stat().st_size, 0)
        custom = self.folder / "explicit_missing_clinvar.pkl.gz"
        self.settings["clinvar_gene_stat"] = str(custom)
        self.write()
        original = self.config.read_bytes()
        with self.assertRaisesRegex(ValueError, "refusing to substitute another dataset/release"):
            self.resolve(output=self.folder / "custom-planning.yaml", source_root=source_root)
        self.assertEqual(self.config.read_bytes(), original)
        with self.assertRaises(ValueError):
            refresh.validate_main(self.settings, source_root, allow_missing_bundled_clinvar=True)
        with mock.patch.object(refresh, "latest_metadata") as network:
            with self.assertRaises((ValueError, OSError)):
                refresh.plan_refresh(self.config, source_root)
        network.assert_not_called()

    def test_explicit_derived_bytes_and_scientific_companion_are_copied_unchanged(self):
        source = self.folder / "selected_motifs.pkl"
        source.write_bytes(b"selected scientific cache bytes")
        sidecar = Path(str(source) + ".nm_to_enst.json")
        sidecar.write_bytes(b'{"mapping_identity":"unchanged fixture"}\n')
        self.settings["alphamissense_intolerant_motifs"] = str(source)
        self.write()
        self.resolve()
        target = Path(yaml.safe_load(self.config.read_text())["alphamissense_intolerant_motifs"])
        self.assertTrue(target.is_relative_to(self.cache))
        self.assertEqual(target.read_bytes(), source.read_bytes())
        self.assertEqual(Path(str(target) + ".nm_to_enst.json").read_bytes(), sidecar.read_bytes())

    def test_preflight_keeps_explicit_source_paths_and_does_not_copy_or_rewrite_selection(self):
        source = self.folder / "exact_selected_reference.fa"
        source.write_bytes(b"selected reference bytes; no sequence validation in path planner")
        self.settings["ref_genome"] = str(source)
        self.write()
        before = self.config.read_bytes()
        output = self.folder / "planning.yaml"
        self.resolve(output=output)
        planning = yaml.safe_load(output.read_text())
        self.assertEqual(planning["ref_genome"], str(source))
        self.assertTrue(Path(planning["alphamissense_pd_stat"]).is_relative_to(self.cache))
        self.assertEqual(self.config.read_bytes(), before)
        self.assertFalse(self.cache.exists())

    def test_checkout_relative_scripts_remain_sources_from_another_working_directory(self):
        self.settings.update(conda_env_yaml="acmg_conda_rough.yml",
                             gene_mechanism_builder_script="scripts/build_gene_pathogenic_mechanism_cache.py")
        self.write()
        previous = Path.cwd()
        try:
            os.chdir(self.folder)
            self.resolve()
        finally:
            os.chdir(previous)
        configured = yaml.safe_load(self.config.read_text())
        for key in ("conda_env_yaml", "gene_mechanism_builder_script"):
            self.assertEqual(configured[key], str(ROOT / self.settings[key]))
        self.assertFalse(self.cache.exists())

    def test_missing_derived_output_is_only_named_under_cache_root(self):
        old_output = self.folder / "old_output/selected_stat.pkl"
        self.settings["alphamissense_pd_stat"] = str(old_output)
        self.write()
        self.resolve()
        target = Path(yaml.safe_load(self.config.read_text())["alphamissense_pd_stat"])
        self.assertTrue(target.is_relative_to(self.cache))
        self.assertFalse(target.exists())
        self.assertFalse(old_output.parent.exists())

    def test_bundle_alias_is_adopted_by_existence_and_size_without_hash_or_copy(self):
        source_root = self.folder / "checkout"
        bundled = source_root / "data/InterPro/Interpro_entry_mapping.pkl.gz"
        bundled.parent.mkdir(parents=True)
        bundled.write_bytes(b"deliberately not a pickle; nonempty bundle policy fixture")
        self.settings["interpro_mapping_pickle"] = bundled.name
        self.write()
        with mock.patch.object(PATHS, "sha256_file", side_effect=AssertionError("bundle content must not be hashed")), \
                mock.patch.object(PATHS.shutil, "copy2", side_effect=AssertionError("bundle copy must not be needed")):
            self.resolve(source_root=source_root)
        self.assertEqual(yaml.safe_load(self.config.read_text())["interpro_mapping_pickle"], str(bundled))
        self.assertFalse(self.cache.exists())

    def test_supplied_source_change_during_copy_does_not_select_partial_copy(self):
        source = self.folder / "selected_scores.vcf.gz"
        source.write_bytes(self.vcf_fixture(body=b"before source bytes"))
        self.settings["spliceai_snv_prescore"] = str(source)
        self.write()
        before = self.config.read_bytes()
        real_copy = PATHS.shutil.copy2
        def mutate(source_path, destination):
            result = real_copy(source_path, destination)
            Path(source_path).write_bytes(b"changed during copy")
            return result
        with mock.patch.object(PATHS.shutil, "copy2", side_effect=mutate):
            with self.assertRaisesRegex(ValueError, "Source changed or copy differs"):
                self.resolve()
        self.assertEqual(self.config.read_bytes(), before)
        self.assertFalse((self.cache / "installer-source-adoptions.json").exists())

    def test_native_spliceai2_tsv_inputs_remain_sources_for_the_manifest_aware_importer(self):
        header = "\t".join(PATHS.spliceai_inspector().native_columns()) + "\n"
        for name in ("renamed_manual_snv.bin", "renamed_manual_indel.vcf.gz"):
            (self.folder / name).write_text(header)
        self.settings.update(spliceai_snv_prescore=str(self.folder / "renamed_manual_snv.bin"),
                             spliceai_indel_prescore=str(self.folder / "renamed_manual_indel.vcf.gz"))
        self.write()
        self.resolve()
        configured = yaml.safe_load(self.config.read_text())
        for key in ("spliceai_snv_prescore", "spliceai_indel_prescore"):
            self.assertEqual(configured[key], self.settings[key])
        self.assertFalse(self.cache.exists())

    def test_completed_spliceai2_pair_uses_the_receipt_aware_physical_relocator(self):
        source = self.folder / "completed_projection"
        source.mkdir()
        names = {role: source / f"SpliceAI2_GRCh38.{role}.vep.vcf.gz" for role in ("snv", "indel")}
        for path in names.values():
            path.write_bytes(self.vcf_fixture(version=2))
            Path(str(path) + ".tbi").write_bytes(b"completed index fixture")
        receipt = source / "adapter.ready.json"
        receipt.write_text(json.dumps({"outputs": {role: {"path": str(path)} for role, path in names.items()}}))
        before = receipt.read_bytes()
        self.settings.update(spliceai_snv_prescore=str(names["snv"]), spliceai_indel_prescore=str(names["indel"]))
        self.write()
        # Resource-specific receipt validation/relocation has its own tests.
        # Here verify delegation, physical ownership and selected config paths.
        def relocate(*, source_dir, cache_dir, cache_root):
            self.assertEqual(source_dir, source)
            self.assertTrue(cache_dir.is_relative_to(cache_root))
            PATHS.shutil.copytree(source_dir, cache_dir)
            (cache_dir / "adapter.original.ready.json").write_bytes(before)
            return {}
        with mock.patch.object(PATHS.spliceai_inspector(), "adopt_completed_projection", side_effect=relocate) as adopter:
            self.resolve()
        adopter.assert_called_once()
        configured = yaml.safe_load(self.config.read_text())
        target = self.cache / "VEP_plugin_cache/SpliceAI2/GRCh38"
        self.assertTrue(target.is_dir())
        self.assertFalse(target.is_symlink())
        for key, role in (("spliceai_snv_prescore", "snv"), ("spliceai_indel_prescore", "indel")):
            self.assertTrue(Path(configured[key]).is_relative_to(self.cache))
            self.assertEqual(Path(configured[key]).resolve(), target / names[role].name)
            self.assertEqual(Path(configured[key]).read_bytes(), names[role].read_bytes())
            self.assertNotEqual(Path(configured[key]).stat().st_ino, names[role].stat().st_ino)
        self.assertEqual((target / "adapter.original.ready.json").read_bytes(), before)
        self.assertEqual(receipt.read_bytes(), before)

    def test_existing_managed_different_bytes_are_retained_before_source_copy(self):
        source = self.folder / "selected_custom.vcf.gz"
        source.write_bytes(self.vcf_fixture(body=b"explicit selected new bytes"))
        self.settings["spliceai_snv_prescore"] = str(source)
        target = Path(PATHS.portable_updates(self.settings, ROOT)["spliceai_snv_prescore"])
        target.parent.mkdir(parents=True)
        target.write_bytes(b"retain old managed bytes")
        self.write()
        self.resolve()
        retained = [path for path in (self.cache / ".installer-retained").rglob("*") if path.is_file()]
        self.assertEqual([path.read_bytes() for path in retained], [b"retain old managed bytes"])
        self.assertEqual(target.read_bytes(), source.read_bytes())

    def test_adoption_never_creates_files_through_a_link_into_an_external_source(self):
        original = self.folder / "unrelated_external_cache"
        original.mkdir()
        self.cache.mkdir()
        (self.cache / "VEP_plugin_cache").symlink_to(original, target_is_directory=True)
        source = self.folder / "selected_input.vcf.gz"
        source.write_bytes(self.vcf_fixture())
        self.settings["spliceai_snv_prescore"] = str(source)
        self.write()
        before = self.config.read_bytes()
        with self.assertRaisesRegex(ValueError, "Refusing to copy through"):
            self.resolve()
        self.assertEqual(list(original.iterdir()), [])
        self.assertEqual(self.config.read_bytes(), before)


if __name__ == "__main__":
    unittest.main()
