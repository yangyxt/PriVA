"""Portable config planning contracts; no downloads or environment mutation."""
import copy
import importlib.util
import json
import os
from pathlib import Path
import re
import subprocess
import sys
import tempfile
import unittest
from unittest import mock

import yaml


ROOT = Path(__file__).resolve().parents[1]
INSTALLER = ROOT / "scripts/install_utils.sh"
HELPER = ROOT / "scripts/installer_config_paths.py"
SPEC = importlib.util.spec_from_file_location("installer_config_paths", HELPER)
PATHS = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(PATHS)
REF_SPEC = importlib.util.spec_from_file_location("portable_reference_resources", ROOT / "scripts/reference_resources.py")
REFERENCES = importlib.util.module_from_spec(REF_SPEC)
REF_SPEC.loader.exec_module(REFERENCES)


class PortableConfigPaths(unittest.TestCase):
    def setUp(self):
        temporary = tempfile.TemporaryDirectory(prefix="priva_portable_config_", dir=os.environ.get("PRIVA_TEST_TMPDIR", os.environ.get("TMPDIR")))
        self.addCleanup(temporary.cleanup)
        self.folder = Path(temporary.name)
        self.cache = self.folder / "selected_cache"
        self.scratch = self.folder / "selected_tmp"
        self.config = self.folder / "selected.yaml"
        self.basic = {"cache_root": str(self.cache), "tmp_dir": str(self.scratch),
                      "assembly": "hg19", "vep_installed_version": "116"}

    def write_config(self, config=None):
        # Keep long test paths on one scalar line so path-syntax cases exercise
        # their intended branch instead of YAML's automatic line wrapping.
        self.config.write_text(yaml.safe_dump(self.basic if config is None else config, sort_keys=False, width=4096))

    def run_helper(self, *extra):
        return subprocess.run([sys.executable, str(HELPER), "resolve", "--config", str(self.config),
                               "--source-root", str(ROOT), *extra], text=True, capture_output=True, timeout=15)

    def run_shell(self, body):
        return subprocess.run(["bash", "--noprofile", "--norc", "-c",
                               'source "$1" >/dev/null || exit\n' + body,
                               "test", str(INSTALLER), str(self.config)],
                              text=True, capture_output=True, timeout=20)

    def test_template_has_no_host_paths_and_only_run_inputs_are_unset(self):
        text = (ROOT / "config.yaml").read_text()
        self.assertNotIn("/paedyl01", text)
        self.assertNotIn("/stornext/", text)
        config = yaml.safe_load(text)
        for key in ("cache_root", "tmp_dir", "input_vcf", "output_dir", "ref_genome"):
            self.assertIn(config[key], PATHS.EMPTY)
        self.assertEqual(str(config["vep_installed_version"]), "116")
        self.assertIsNone(config["nm_to_enst_tsv"])
        self.assertIn("not for installation", text)

    def test_blank_basic_paths_stop_before_manager_or_any_config_mutation(self):
        self.write_config(dict(self.basic, cache_root="", tmp_dir=""))
        before = self.config.read_bytes()
        result = self.run_shell('''
            conda_install_vep() { printf 'UNEXPECTED_MANAGER\n'; return 99; }
            main_install "$2"
        ''')
        self.assertNotEqual(result.returncode, 0)
        self.assertIn("set cache_root", result.stderr)
        self.assertIn("set tmp_dir", result.stderr)
        self.assertNotIn("UNEXPECTED", result.stdout)
        self.assertEqual(self.config.read_bytes(), before)
        self.assertFalse(self.scratch.exists())

    def test_required_basic_gate_rejects_invalid_storage_assembly_and_release(self):
        self.cache.mkdir()
        self.scratch.mkdir()
        for change, diagnostic in (
                ({"cache_root": "/"}, "set cache_root to an absolute, dedicated storage directory"),
                ({"tmp_dir": "relative/tmp"}, "set tmp_dir to an absolute, dedicated storage directory"),
                ({"cache_root": str(self.folder / "unsafe space")}, "without spaces or shell characters"),
                ({"assembly": "unknown"}, "set assembly"),
                ({"vep_installed_version": "116.3"}, "set vep_installed_version")):
            with self.subTest(change=change):
                self.write_config(dict(self.basic, **change))
                before = self.config.read_bytes()
                result = self.run_shell('priva_install_basic_paths --config "$2"')
                self.assertNotEqual(result.returncode, 0)
                self.assertIn(diagnostic, result.stderr)
                self.assertEqual(self.config.read_bytes(), before)
                self.assertEqual(list(self.cache.iterdir()), [])
                self.assertEqual(list(self.scratch.iterdir()), [])

    def test_missing_or_nondirectory_storage_stops_before_environment_or_config_changes(self):
        for state in ("missing", "regular_file", "broken_symlink"):
            for invalid_keys in (("cache_root",), ("tmp_dir",), ("cache_root", "tmp_dir")):
                for entry_point in ("priva_install_basic_paths --config", "main_install"):
                    with self.subTest(state=state, invalid_keys=invalid_keys, entry_point=entry_point):
                        case = self.folder / (state + "_" + "_".join(invalid_keys) + "_" + entry_point.split()[0])
                        case.mkdir()
                        paths = {key: case / key / "selected" for key in ("cache_root", "tmp_dir")}
                        for key, path in paths.items():
                            if key not in invalid_keys:
                                path.mkdir(parents=True)
                            elif state == "regular_file":
                                path.parent.mkdir()
                                path.write_bytes(b"selected path is an existing file\n")
                            elif state == "broken_symlink":
                                path.parent.mkdir()
                                path.symlink_to(path.parent / "missing_target", target_is_directory=True)
                        self.write_config(dict(self.basic, **{key: str(path) for key, path in paths.items()}))
                        before = self.config.read_bytes()
                        entries = sorted(str(path.relative_to(case)) for path in case.rglob("*"))
                        result = self.run_shell('''
                            conda_install_vep() { printf 'UNEXPECTED_MANAGER\n'; return 91; }
                            priva_prepare_install_config() { printf 'UNEXPECTED_CONFIG_PREPARATION\n'; return 92; }
                            priva_configure_cache_paths() { printf 'UNEXPECTED_CACHE_PREPARATION\n'; return 93; }
                            priva_configure_portable_paths() { printf 'UNEXPECTED_CONFIG_WRITE\n'; return 94; }
                            reference_genome_install() { printf 'UNEXPECTED_RESOURCE\n'; return 95; }
                            wget() { printf 'UNEXPECTED_DOWNLOAD\n'; return 96; }
                            curl() { printf 'UNEXPECTED_DOWNLOAD\n'; return 96; }
                        ''' + '\n' + entry_point + ' "$2"')
                        self.assertNotEqual(result.returncode, 0)
                        for key in paths:
                            diagnostic = f"{key} must be an existing directory"
                            if key in invalid_keys:
                                self.assertIn(diagnostic, result.stderr)
                            else:
                                self.assertNotIn(diagnostic, result.stderr)
                        self.assertNotIn("UNEXPECTED", result.stdout + result.stderr)
                        self.assertEqual(self.config.read_bytes(), before)
                        self.assertEqual(sorted(str(path.relative_to(case)) for path in case.rglob("*")), entries)
                        for key in invalid_keys:
                            path = paths[key]
                            if state == "regular_file":
                                self.assertEqual(path.read_bytes(), b"selected path is an existing file\n")
                            elif state == "broken_symlink":
                                self.assertTrue(path.is_symlink())
                                self.assertEqual(path.readlink(), path.parent / "missing_target")
                                self.assertFalse(path.exists())
                            else:
                                self.assertFalse(path.parent.exists())

    def test_basic_gate_needs_no_yaml_environment_and_preserves_selected_temp(self):
        self.cache.mkdir()
        self.scratch.mkdir()
        self.write_config()
        before = self.config.read_bytes()
        result = self.run_shell('''
            python() { printf 'UNEXPECTED_PYTHON\n'; return 99; }
            yq() { printf 'UNEXPECTED_YQ\n'; return 99; }
            priva_install_basic_paths --config "$2" || exit
            printf '%s\n' "$PRIVA_INSTALL_TMPDIR"
        ''')
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertEqual(result.stdout.strip(), str(self.scratch))
        self.assertEqual(self.config.read_bytes(), before)
        self.assertEqual(list(self.cache.iterdir()), [])
        self.assertEqual(list(self.scratch.iterdir()), [])

    def test_existing_directory_symlinks_pass_basic_gate_without_rewriting_paths(self):
        actual_cache = self.folder / "actual_cache"
        actual_tmp = self.folder / "actual_tmp"
        actual_cache.mkdir()
        actual_tmp.mkdir()
        self.cache.symlink_to(actual_cache, target_is_directory=True)
        self.scratch.symlink_to(actual_tmp, target_is_directory=True)
        self.write_config()
        before = self.config.read_bytes()
        result = self.run_shell('''
            priva_install_basic_paths --config "$2" || exit
            printf '%s\n' "$PRIVA_INSTALL_TMPDIR"
        ''')
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertEqual(result.stdout.strip(), str(self.scratch))
        self.assertEqual(self.config.read_bytes(), before)
        self.assertTrue(self.cache.is_symlink())
        self.assertTrue(self.scratch.is_symlink())
        self.assertEqual(list(actual_cache.iterdir()), [])
        self.assertEqual(list(actual_tmp.iterdir()), [])

    def test_main_reaches_environment_with_existing_storage_without_mkdir(self):
        self.cache.mkdir()
        self.scratch.mkdir()
        self.write_config()
        before = self.config.read_bytes()
        result = self.run_shell('''
            mkdir() { printf 'UNEXPECTED_DIRECTORY_CREATION\n'; return 99; }
            conda_install_vep() { printf 'ENVIRONMENT_REACHED %s\n' "$TMPDIR"; return 73; }
            priva_prepare_install_config() { printf 'UNEXPECTED_CONFIG_PREPARATION\n'; return 99; }
            main_install "$2"
        ''')
        self.assertNotEqual(result.returncode, 0)
        self.assertIn("ENVIRONMENT_REACHED " + str(self.scratch), result.stdout)
        self.assertNotIn("UNEXPECTED", result.stdout + result.stderr)
        self.assertEqual(self.config.read_bytes(), before)
        self.assertEqual(list(self.cache.iterdir()), [])
        self.assertEqual(list(self.scratch.iterdir()), [])

    def test_minimal_config_plans_both_assemblies_and_aliases(self):
        for assembly, ucsc, genome in (("hg19", "hg19", "GRCh37"), ("GRCh37", "hg19", "GRCh37"),
                                       ("hg38", "hg38", "GRCh38"), ("GRCh38", "hg38", "GRCh38")):
            with self.subTest(assembly=assembly):
                config = dict(self.basic, assembly=assembly)
                updates = PATHS.portable_updates(config, ROOT)
                self.assertEqual(updates["clinvar_vcf"], str(self.cache / f"ClinVar/clinvar.{ucsc}.vep.vcf.gz"))
                self.assertEqual(updates["cds_fasta_file"], str(self.cache / f"VEP_plugin_cache/Ensembl_CDS/release-116/Homo_sapiens.{genome}.cds.all.fa.gz"))
                self.assertEqual(updates["repeat_region_file_name"], f"ucsc.{ucsc}.fasta.repeat.bed.gz")
                self.assertEqual(updates["loftee_repo"], str(self.cache / f"VEP_plugin/LoFtee_repos/loftee-{assembly}"))
                resources, _ = REFERENCES.resource_plan(dict(config, **updates), ROOT / "data", "gofcards", use_cache_root=True)
                paths = {str(resource["path"]) for resource in resources}
                self.assertIn(str(self.cache / "reference/ucsc.hg19.fasta"), paths)
                self.assertIn(str(self.cache / "reference/ucsc.hg38.fasta"), paths)
                self.assertIn(str(self.cache / "liftOver/hg19ToHg38.over.chain"), paths)
                self.assertFalse(self.cache.exists())

    def test_generated_destinations_use_storage_and_bundled_sources_use_checkout(self):
        updates = PATHS.portable_updates(self.basic, ROOT)
        for key in ("gene_mechanism_cache_dir", "gene_mechanism_raw_dir",
                    "gofcards_workdir", "clinvar_vcv_xml", "hpo_raw_dir", "hpo_assertions",
                    "mondo_obo", "mondo_disease_scope_registry"):
            with self.subTest(key=key):
                self.assertTrue(Path(updates[key]).is_relative_to(self.cache))
                self.assertFalse(Path(updates[key]).is_relative_to(ROOT))
        for key, relative in PATHS.BUNDLED_PATHS.items():
            self.assertEqual(updates[key], str(ROOT / relative))
        self.assertEqual(updates["interpro_mapping_pickle"], str(ROOT / "data/InterPro/Interpro_entry_mapping.pkl.gz"))
        for key, relative in PATHS.BUNDLED_CACHES.items():
            self.assertEqual(updates[key], str(ROOT / relative))
        self.assertEqual(updates["clingen_map"], str(ROOT / "data/clingen/clingen_map.hg19.pkl.gz"))
        self.assertEqual(updates["gene_dosage_sensitivity"], str(ROOT / "data/clingen/gene_dosage_sensitivity.hg19.tsv"))
        for key in ("ref_genome", "gofcards_hg19_fasta", "gofcards_hg38_fasta", "gofcards_hg19_to_hg38_chain",
                    "nm_to_enst_tsv", "primateai_prescore", "spliceai_snv_prescore", "spliceai_indel_prescore",
                    "input_vcf", "output_dir", "mavedb_file", "pext_mean_bw"):
            self.assertNotIn(key, updates)

    def test_external_cache_sources_get_managed_destinations_while_scripts_remain_sources(self):
        external = self.folder / "selected external's input"
        config = dict(self.basic, interpro_mapping_pickle="Interpro_entry_mapping.pkl.gz",
                      clingen_map=str(external), hpo_assertions=str(external), ref_genome=str(external),
                      cds_fasta_file=str(external), gene_mechanism_cache_dir=str(external),
                      gene_mechanism_builder_script=str(external), nm_to_enst_tsv=str(external))
        updates = PATHS.portable_updates(config, ROOT)
        for key in ("clingen_map", "hpo_assertions", "ref_genome", "cds_fasta_file", "gene_mechanism_cache_dir", "nm_to_enst_tsv"):
            self.assertTrue(Path(updates[key]).is_relative_to(self.cache))
        self.assertNotIn("gene_mechanism_builder_script", updates)
        self.assertTrue(Path(updates["gene_mechanism_raw_dir"]).is_relative_to(self.cache))
        self.assertFalse(external.exists())

    def test_existing_managed_root_migration_is_applied_before_defaults(self):
        source = INSTALLER.read_text().split("<<'PRIVA_CACHE_PY'\n", 1)[1].split("\nPRIVA_CACHE_PY", 1)[0]
        namespace = {"__name__": "portable_test_managed_resolver", "CONFIG_PATHS_HELPER": ROOT / "scripts/installer_config_paths.py"}
        exec(compile(source, str(INSTALLER), "exec"), namespace)
        old = self.folder / "old_cache"
        external = self.folder / "external.vcf.gz"
        config = dict(self.basic, vep_plugins_cachedir=str(old / "plugins"),
                      human_ancestor_fasta=str(old / "plugins/LoFtee/hg19/human_ancestor.fa.gz"),
                      clinvar_vcf=str(external), cache_folders={"vep_plugin_cache": "plugins_new"})
        _, managed = namespace["resolve_cache_paths"](config)
        effective = dict(config, **managed)
        updates = PATHS.portable_updates(effective, ROOT)
        effective.update(updates)
        self.assertEqual(effective["human_ancestor_fasta"], str(self.cache / "plugins_new/LoFtee/hg19/human_ancestor.fa.gz"))
        self.assertTrue(Path(effective["clinvar_vcf"]).is_relative_to(self.cache))
        self.assertEqual(effective["splicevault_prescore"], str(self.cache / "plugins_new/SpliceVault/SpliceVault_data_GRCh37.tsv.gz"))
        self.assertFalse(old.exists())
        self.assertFalse(self.cache.exists())

    def test_relative_bundled_paths_resolve_independently_of_working_directory(self):
        config = dict(self.basic, hgnc_table="data/hgnc/non_alt_loci_set.tsv",
                      gene_mechanism_builder_script="scripts/build_gene_pathogenic_mechanism_cache.py")
        before = copy.deepcopy(config)
        updates = PATHS.portable_updates(config, ROOT)
        self.assertEqual(updates["hgnc_table"], str(ROOT / config["hgnc_table"]))
        self.assertEqual(updates["gene_mechanism_builder_script"], str(ROOT / config["gene_mechanism_builder_script"]))
        self.assertEqual(config, before)

    def test_check_output_is_read_only_and_write_is_idempotent_with_comments(self):
        self.write_config()
        self.config.write_text("# selected config\n" + self.config.read_text() + "hpo_assertions: ''  # keep this comment\nthreads: 20\n")
        before = self.config.read_bytes()
        result = self.run_helper("--check")
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertEqual(json.loads(result.stdout)["hpo_assertions"], str(self.cache / "hpo/genes_to_phenotype.assertions.tsv.gz"))
        self.assertEqual(self.config.read_bytes(), before)
        self.assertFalse(self.cache.exists())
        self.assertEqual(self.run_helper().returncode, 0)
        after = self.config.read_bytes()
        self.assertIn(b"# selected config", after)
        self.assertIn(b"# keep this comment", after)
        self.assertIn(b"threads: 20", after)
        self.assertEqual(self.run_helper().returncode, 0)
        self.assertEqual(self.config.read_bytes(), after)
        self.assertFalse(self.cache.exists())

    def test_effective_preflight_output_keeps_source_config_and_applies_managed_updates(self):
        self.write_config()
        before = self.config.read_bytes()
        output = self.folder / "effective.yaml"
        selected = str(self.folder / "other_selected_plugins")
        result = self.run_helper("--managed-updates-json", json.dumps({"vep_plugins_cachedir": selected}), "--output", str(output))
        self.assertNotEqual(result.returncode, 0)
        self.assertIn("differ", result.stderr)
        self.assertEqual(self.config.read_bytes(), before)
        self.assertFalse(output.exists())
        self.assertFalse(self.cache.exists())

    def test_missing_custom_reference_is_not_replaced_by_managed_default(self):
        missing = self.folder / "selected_missing_reference.fa"
        config = dict(self.basic, ref_genome=str(missing))
        updates = PATHS.portable_updates(config, ROOT)
        with self.assertRaisesRegex(ValueError, "Explicit.*missing"):
            PATHS.adoption_plan(config, updates, ROOT)
        self.assertEqual(config["ref_genome"], str(missing))
        self.assertFalse(self.cache.exists())

    def test_duplicate_keys_fail_before_writing_resolved_config(self):
        self.write_config()
        self.config.write_text(self.config.read_text() + 'hpo_assertions: one\nhpo_assertions: two\n')
        before = self.config.read_bytes()
        result = self.run_helper()
        self.assertNotEqual(result.returncode, 0)
        self.assertIn("Duplicate configuration key", result.stderr)
        self.assertEqual(self.config.read_bytes(), before)

    def test_missing_or_external_interpro_bundle_keeps_downloads_in_storage(self):
        source = self.folder / "checkout"
        bundled = source / "data/InterPro/Interpro_entry_mapping.pkl.gz"
        self.assertEqual(PATHS.portable_updates(self.basic, source)["interpro_mapping_pickle"],
                         str(self.cache / "InterPro/Interpro_entry_mapping.pkl.gz"))
        bundled.parent.mkdir(parents=True)
        external = self.folder / "external_interpro.pkl.gz"
        external.write_bytes(b"supplied bytes")
        bundled.symlink_to(external)
        self.assertEqual(PATHS.portable_updates(self.basic, source)["interpro_mapping_pickle"],
                         str(self.cache / "InterPro/Interpro_entry_mapping.pkl.gz"))
        self.assertEqual(external.read_bytes(), b"supplied bytes")

    def test_bundled_interpro_needs_only_nonzero_file_without_xml_or_decode(self):
        source = self.folder / "checkout"
        bundled = source / "data/InterPro/Interpro_entry_mapping.pkl.gz"
        bundled.parent.mkdir(parents=True)
        bundled.write_bytes(b"existence policy fixture, deliberately not a pickle")
        self.write_config()
        result = self.run_shell('''
            DATA_DIR="$(dirname "$2")/checkout/data"
            read_yaml() { [[ $2 != interpro_mapping_pickle ]] || printf 'Interpro_entry_mapping.pkl.gz\n'; }
            update_yaml() { printf 'REGISTER %s\n' "$3"; }
            update_or_append_yaml() { update_yaml "$@"; }
            wget() { printf 'UNEXPECTED_DOWNLOAD\n'; return 99; }
            gzip() { printf 'UNEXPECTED_CONTENT_CHECK\n'; return 99; }
            python() { printf 'UNEXPECTED_DECODE\n'; return 99; }
            InterPro_parsing "$2"
        ''')
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertIn("REGISTER " + str(bundled), result.stdout)
        self.assertNotIn("UNEXPECTED", result.stdout)
        self.assertEqual(bundled.read_bytes(), b"existence policy fixture, deliberately not a pickle")
        self.assertEqual(sorted(p.name for p in bundled.parent.iterdir()), [bundled.name])

    def test_known_bundled_clinvar_uses_size_policy_but_external_symlink_stays_strict(self):
        spec = importlib.util.spec_from_file_location("portable_bundle_prerequisites", ROOT / "scripts/installer_prerequisites.py")
        prereqs = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(prereqs)
        source = self.folder / "checkout"
        bundled = source / "data/ClinVar/clinvar_2star_stats.pkl.gz"
        bundled.parent.mkdir(parents=True)
        bundled.write_bytes(b"nonempty original bundle fixture")
        with mock.patch.object(prereqs, "validate_clinvar_stats", side_effect=AssertionError("content check")):
            prereqs.validate_selected_clinvar_stats({}, source)
        bundled.unlink()
        external = self.folder / "external_stats.pkl.gz"
        external.write_bytes(b"nonempty external invalid pickle")
        bundled.symlink_to(external)
        with self.assertRaisesRegex(ValueError, "Cannot reuse bundled ClinVar"):
            prereqs.validate_selected_clinvar_stats({}, source)

    def test_main_uses_effective_destinations_but_keeps_source_config_on_custom_failure(self):
        self.cache.mkdir()
        self.scratch.mkdir()
        self.write_config(dict(self.basic, hpo_assertions=str(self.folder / "explicit.hpo.gz")))
        before = self.config.read_bytes()
        result = self.run_shell('''
            conda_install_vep() { PRIVA_INSTALL_ENV_NAME=priva_acmg; }
            # Labelled boundary stub: this fixture tests the effective config
            # and custom prerequisite failure; real manual-input gating has
            # separate fresh/recovery integration coverage.
            priva_manual_cache_inputs() { return 0; }
            priva_configure_cache_paths() {
                [[ ${2:-} == --check ]] || { printf 'UNEXPECTED_CACHE_WRITE\n'; return 97; }
                printf '{}\n'
            }
            priva_install_prerequisites() {
                python - "$1" <<'PY'
import sys,yaml
from pathlib import Path
c=yaml.safe_load(Path(sys.argv[1]).read_text())
print('EFFECTIVE',c['gofcards_workdir'],c['hpo_assertions'],c['gofcards_builder_script'])
PY
                return 73
            }
            nm_to_enst_mapping_install() { printf 'UNEXPECTED_MAPPING\n'; return 99; }
            main_install "$2"
        ''')
        self.assertNotEqual(result.returncode, 0)
        self.assertIn("EFFECTIVE " + str(self.cache / "gofcards/build_work"), result.stdout)
        self.assertIn(str(self.cache / "hpo/supplied/hpo_assertions/explicit.hpo.gz"), result.stdout)
        self.assertIn(str(ROOT / "scripts/build_gofcards_exact_gof_cache.py"), result.stdout)
        self.assertNotIn("UNEXPECTED", result.stdout)
        self.assertEqual(self.config.read_bytes(), before)
        self.assertFalse(list(self.scratch.iterdir()))


if __name__ == "__main__":
    unittest.main()
