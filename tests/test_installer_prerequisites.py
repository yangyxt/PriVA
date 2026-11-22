"""Original-input and call contracts; fixture-only, without resource downloads."""
import gzip
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


ROOT = Path(__file__).resolve().parents[1]
INSTALLER = ROOT / "scripts/install_utils.sh"
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
                        PRIVA_INSTALL_INSPECT_ROOTS=os.pathsep.join((str(ROOT), str(self.folder))))
        self.env.pop("PRIVA_SMOKE_SKIP_VEP", None)

    def write_config(self):
        self.config_file.write_text(json.dumps(self.config))

    def run_shell(self, body):
        self.write_config()
        # Avoid a yq dependency in these shell fixtures while retaining the
        # actual Python prerequisite checks and command argument boundaries.
        prefix = '''
            source "$1"
            log() { printf '%s\\n' "$*" >&2; }
            # This module tests prerequisites; basic YAML bootstrap is separate.
            priva_install_basic_paths() { PRIVA_INSTALL_TMPDIR="$TEST_FOLDER"; }
            # Labelled main-install boundary stub: real manual-gate and prompt
            # checks live in test_install_fresh_recovery_integration.py.
            priva_manual_cache_inputs() { return 0; }
            read_yaml() {
                "$TEST_PYTHON" -c 'import json,sys; v=json.load(open(sys.argv[1])).get(sys.argv[2]); print("null" if v is None else v)' "$1" "$2"
            }
            update_yaml() { printf 'UPDATE:%s\\n' "$2"; }
            update_or_append_yaml() { update_yaml "$@"; }
            wget() { printf 'UNEXPECTED_DOWNLOAD\\n'; return 97; }
            curl() { printf 'UNEXPECTED_DOWNLOAD\\n'; return 97; }
            fixture_annotation_base() {
                # This module isolates prerequisite/mapping routing. Complete
                # base science/publication has its own real BGZF fixtures.
                case $2 in
                    path) "$TEST_PYTHON" -c 'import json,sys;print(json.load(open(sys.argv[1]))["alphamissense_vep_vcf"])' "$TEST_CONFIG" ;;
                    check) return 0 ;;
                    motif-check) [[ -s ${@: -1} ]] ;;
                    *) return 99 ;;
                esac
            }
            python() {
                if [[ $1 == */alphamissense_annotation_base.py ]]; then
                    fixture_annotation_base "$@"
                else "$TEST_PYTHON" "$@"; fi
            }
        '''
        if "python() {" in body:
            # Keep each labelled scientific stub across the publisher's child
            # process boundary; exercise the real transaction helper itself.
            body = body.replace("python() {", "fixture_python() {")
            self.env["TEST_PUBLISHER_DRIVER"] = str(ROOT / "tests/mapped_artifact_fixture_driver.py")
            prefix += r'''
                python() {
                    if [[ $1 == */alphamissense_annotation_base.py ]]; then
                        fixture_annotation_base "$@"
                    elif [[ $1 == */mapped_artifact_publish.py ]]; then
                        export -f fixture_python
                        "$TEST_PYTHON" "$TEST_PUBLISHER_DRIVER" run --publisher "$1" -- "${@:2}"
                    else fixture_python "$@"; fi
                }
            '''
        return subprocess.run(["bash", "--noprofile", "--norc", "-c", prefix + body,
                               "test", str(INSTALLER)], env=self.env, text=True,
                              capture_output=True, timeout=30)

    def test_explicit_missing_mapping_stops_main_before_downloads_or_config_writes(self):
        self.config.update(cache_root=str(self.folder / "new_cache"), tmp_dir=str(self.folder))
        self.config["nm_to_enst_tsv"] = str(self.folder / "missing_original_mapping.tsv")
        result = self.run_shell('''
            conda_install_vep() { PRIVA_INSTALL_ENV_NAME=priva_acmg; }
            priva_configure_cache_paths() { [[ ${2:-} != --check ]] || { printf '{}\\n'; return; }; printf 'UNEXPECTED_CONFIG_WRITE\\n'; }
            reference_genome_install() { printf 'UNEXPECTED_REFERENCE\\n'; }
            main_install "$TEST_CONFIG"
        ''')
        self.assertNotEqual(result.returncode, 0)
        self.assertIn("nm_to_enst_tsv", result.stderr)
        self.assertNotIn("UNEXPECTED", result.stdout)
        self.assertEqual(json.loads(self.config_file.read_text()), self.config)
        self.assertFalse((self.folder / "build_work").exists())

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

    def test_clinvar_installer_registers_nonempty_original_bundle_without_content_validation(self):
        self.stats.write_bytes(b"nonempty original bundle fixture")
        before = self.stats.read_bytes()
        result = self.run_shell('''
            BASE_DIR="$TEST_FOLDER/selected_source"
            DATA_DIR="$BASE_DIR/data"
            update_or_append_yaml() { printf 'DEPLOY:%s\\n' "$3"; }
            ClinVar_Gene_stat "$TEST_CONFIG"
        ''')
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertEqual(result.stdout.splitlines(), ["DEPLOY:" + str(self.stats)])
        self.assertEqual(self.stats.read_bytes(), before)
        self.stats.write_bytes(b"")
        result = self.run_shell('''
            BASE_DIR="$TEST_FOLDER/selected_source"
            DATA_DIR="$BASE_DIR/data"
            python() {
                [[ $1 == */installer_prerequisites.py ]] || { printf 'UNEXPECTED_SCIENTIFIC_BUILD\\n'; return 99; }
                "$TEST_PYTHON" "$@"
            }
            ClinVar_Gene_stat "$TEST_CONFIG"
        ''')
        self.assertNotEqual(result.returncode, 0)
        self.assertIn("ClinVar_Gene_stat_regenerate CONFIG SUMMARY MANE OUTPUT", result.stderr)
        self.assertNotIn("UNEXPECTED", result.stdout)
        self.assertNotIn("UPDATE", result.stdout)

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

    def test_main_preflight_preserves_requested_release_without_config_mutation(self):
        self.config.update(cache_root=str(self.folder / "new_cache"), tmp_dir=str(self.folder))
        Path(self.config["cache_root"]).mkdir()
        self.config["vep_installed_version"] = "999"
        for assembly in ("GRCh37", "GRCh38"):
            shutil.copytree(self.cache / ("homo_sapiens_merged/113_" + assembly),
                            self.cache / ("homo_sapiens_merged/999_" + assembly))
        result = self.run_shell('''
            conda_install_vep() { PRIVA_INSTALL_ENV_NAME=priva_acmg; }
            priva_verify_vep_release() { [[ $1 == 999 ]]; }
            priva_configure_cache_paths() { [[ ${2:-} != --check ]] || { printf '{}\\n'; return; }; printf 'CONFIG_STAGE_REACHED\\n'; return 73; }
            reference_genome_install() { printf 'UNEXPECTED_REFERENCE\\n'; }
            main_install "$TEST_CONFIG"
        ''')
        self.assertIn("CONFIG_STAGE_REACHED", result.stdout, result.stderr)
        self.assertNotIn("UNEXPECTED", result.stdout)
        self.assertEqual(json.loads(self.config_file.read_text())["vep_installed_version"], "999")
        self.assertFalse((self.folder / "build_work").exists())

    def test_mechanism_source_injection_and_build_order_is_preserved(self):
        result = self.run_shell('''
            mondo_hpo_scope_install() { printf 'MONDO\\n'; }
            gene_pathogenic_mechanism_cache_install() { printf 'GENE\\n'; }
            gofcards_exact_gof_cache_install() { printf 'NORMALIZE\\n'; }
            gene_nonlof_mechanism_cache_install() { printf 'NONLOF:%s\\n' "${2:-build}"; }
            gofcards_clinvar_injection_install() { printf 'INJECT\\n'; }
            hpo_condition_mechanism_cache_install() { printf 'HPO\\n'; }
            mechanism_resource_install "$TEST_CONFIG"
        ''')
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertEqual(result.stdout.splitlines(), [
            "MONDO", "GENE", "NORMALIZE", "NONLOF:sources-only", "INJECT", "NONLOF:build", "HPO",
        ])

    def test_mapping_is_passed_verbatim_and_pm1_options_are_preserved(self):
        self.config.update(cache_root=str(self.folder), vep_plugins_cachedir=str(self.folder))
        result = self.run_shell('''
            python() {
                if [[ $1 == */installer_prerequisites.py || $1 == */biomart_nm_to_enst.py ]]; then
                    "$TEST_PYTHON" "$@"
                else
                    "$TEST_PYTHON" -c 'import json,sys; from pathlib import Path; args=sys.argv[1:]; print("CALL:"+json.dumps(args)); Path(args[args.index("--output")+1]).write_text("fixture scientific output\\n")' "$@"
                fi
            }
            AlphaMissense_pick_missense_constrained_segments "$TEST_CONFIG" || exit
            prepare_pm1_regions "$TEST_CONFIG"
        ''')
        self.assertEqual(result.returncode, 0, result.stderr)
        calls = [json.loads(line[5:]) for line in result.stdout.splitlines() if line.startswith("CALL:")]
        self.assertEqual(len(calls), 2)
        for call in calls:
            self.assertEqual(call[call.index("--nm_to_enst_tsv") + 1], str(self.mapping))
        self.assertEqual(calls[1][calls[1].index("--threads") + 1], "2")
        self.assertIn("--skip_das", calls[1])
        self.assertNotIn("--skip_das", calls[0])

    def test_assembly_aliases_match_canonical_motif_and_pm1_calls(self):
        self.config.update(cache_root=str(self.folder), vep_plugins_cachedir=str(self.folder))
        body = r'''
            python() {
                if [[ $1 == */installer_prerequisites.py || $1 == */biomart_nm_to_enst.py ]]; then
                    "$TEST_PYTHON" "$@"
                else
                    "$TEST_PYTHON" -c '
import json,sys
from pathlib import Path
args=sys.argv[1:]
for option in ("--rmc_tsv", "--hcseeker_tsv"):
    if option in args:
        assert Path(args[args.index(option)+1]).is_file()
if "--assembly" in args:
    assert args[args.index("--assembly")+1] in ("hg19", "hg38")
print("CALL:"+json.dumps(args))
Path(args[args.index("--output")+1]).write_bytes(b"fixture scientific output")
' "$@"
                fi
            }
            AlphaMissense_pick_missense_constrained_segments "$TEST_CONFIG" || exit
            prepare_pm1_regions "$TEST_CONFIG"
        '''
        for canonical, alias in (("hg19", "GRCh37"), ("hg38", "GRCh38")):
            observed = []
            for assembly in (canonical, alias):
                with self.subTest(assembly=assembly):
                    for path in (self.folder / "motifs.pkl", self.folder / f"pm1_regions.{canonical}.pkl"):
                        path.unlink(missing_ok=True)
                        Path(str(path) + ".nm_to_enst.json").unlink(missing_ok=True)
                    self.config["assembly"] = assembly
                    result = self.run_shell(body)
                    self.assertEqual(result.returncode, 0, result.stderr)
                    calls = [json.loads(line[5:]) for line in result.stdout.splitlines() if line.startswith("CALL:")]
                    self.assertEqual(len(calls), 2)
                    self.assertEqual(calls[0][calls[0].index("--rmc_tsv") + 1],
                                     str(ROOT / f"scripts/../data/gnomAD_RMC/gnomad_v2.1.1_rmc_flat.{canonical}.tsv"))
                    self.assertEqual(calls[0][calls[0].index("--hcseeker_tsv") + 1],
                                     str(ROOT / f"scripts/../data/HCSeeker/HC_spots.{canonical}.tsv"))
                    self.assertEqual(calls[1][calls[1].index("--assembly") + 1], canonical)
                    for call in calls:
                        index = call.index("--output") + 1
                        staged = Path(call[index])
                        self.assertTrue(staged.parent.name.startswith(".mapped-artifact-"))
                        final = self.folder / staged.name
                        self.assertTrue(final.is_file())
                        call[index] = str(final)
                    self.assertEqual(calls[1][calls[1].index("--output") + 1],
                                     str(self.folder / f"pm1_regions.{canonical}.pkl"))
                    self.assertIn("--skip_das", calls[1])
                    self.assertEqual(json.loads(self.config_file.read_text())["assembly"], assembly)
                    observed.append(calls)
            self.assertEqual(observed[0], observed[1])

    def test_assembly_aliases_match_canonical_domain_filters_and_destinations(self):
        Path(self.config["alphamissense_vep_vcf"]).write_bytes(b"fixture VCF")
        self.config["clinvar_vcf"] = str(self.folder / "clinvar.vcf.gz")
        Path(self.config["clinvar_vcf"]).write_bytes(b"fixture ClinVar")
        Path(self.config["alphamissense_pd_stat"]).write_bytes(b"opaque stats fixture")
        self.config["interpro_mapping_pickle"] = str(self.folder / "mapping.pkl")
        Path(self.config["interpro_mapping_pickle"]).write_bytes(b"opaque mapping fixture")
        body = r'''
            rm -f "$TEST_FOLDER"/domain_tolerance_analysis.* "$TEST_FOLDER"/transcript_pp_domain_mapping.* "$TEST_FOLDER"/all_intolerant_domains.*
            python() {
                if [[ $1 != */am_pick_intolerant_domains.py ]]; then "$TEST_PYTHON" "$@"; return $?; fi
                "$TEST_PYTHON" -c '
import json,sys
from pathlib import Path
args=sys.argv[1:]
folder=Path(args[args.index("--output_dir")+1])
assembly=args[args.index("--assembly")+1]
for stem,suffix in (("domain_tolerance_analysis","tsv"),("transcript_pp_domain_mapping","pkl"),("all_intolerant_domains","pkl")):
    (folder / f"{stem}.{assembly}.{suffix}").write_bytes(b"fixture scientific output")
print("CALL:"+json.dumps(args))
' "$@"; }
            update_yaml() { printf 'DEPLOY:%s=%s\n' "$2" "$3"; }
            AlphaMissense_pick_intolerant_domains "$TEST_CONFIG"
        '''
        for canonical, alias in (("hg19", "GRCh37"), ("hg38", "GRCh38")):
            observed = []
            for assembly in (canonical, alias):
                with self.subTest(assembly=assembly):
                    self.config["assembly"] = assembly
                    result = self.run_shell(body)
                    self.assertEqual(result.returncode, 0, result.stderr)
                    calls = [json.loads(line[5:]) for line in result.stdout.splitlines() if line.startswith("CALL:")]
                    self.assertEqual(len(calls), 1)
                    call = calls[0]
                    self.assertEqual(call[call.index("--assembly") + 1], canonical)
                    self.assertEqual(call[call.index("--hcseeker_spots_tsv") + 1],
                                     str(ROOT / f"data/HCSeeker/HC_spots.{canonical}.tsv"))
                    self.assertEqual(call[call.index("--am_vcf") + 1], self.config["alphamissense_vep_vcf"])
                    self.assertEqual(call[call.index("--clinvar_vcf") + 1], self.config["clinvar_vcf"])
                    self.assertIn(f"domain_tolerance_analysis.{canonical}.tsv", result.stdout)
                    self.assertIn(f"transcript_pp_domain_mapping.{canonical}.pkl", result.stdout)
                    self.assertIn(f"all_intolerant_domains.{canonical}.pkl", result.stdout)
                    self.assertEqual(json.loads(self.config_file.read_text())["assembly"], assembly)
                    observed.append(result.stdout)
            self.assertEqual(observed[0], observed[1])

    def test_existing_motif_or_pm1_without_mapping_receipt_is_preserved(self):
        for function, name in (("AlphaMissense_pick_missense_constrained_segments", "motifs.pkl"),
                               ("prepare_pm1_regions", "pm1_regions.hg38.pkl")):
            with self.subTest(function=function):
                artifact = self.folder / name
                artifact.write_bytes(b"existing artifact with unknown mapping\n")
                before = artifact.read_bytes()
                result = self.run_shell(f'{function} "$TEST_CONFIG"')
                self.assertNotEqual(result.returncode, 0)
                self.assertIn("missing/stale NM-to-ENST receipt", result.stderr)
                self.assertEqual(artifact.read_bytes(), before)
                self.assertNotIn("UPDATE", result.stdout)

    def test_existing_fallback_motif_cannot_bypass_receipt_guard(self):
        self.config["alphamissense_intolerant_motifs"] = None
        artifact = self.folder / "fallback.pkl"
        artifact.write_bytes(b"existing fallback artifact\n")
        result = self.run_shell('AlphaMissense_pick_missense_constrained_segments "$TEST_CONFIG" "$TEST_FOLDER/fallback"')
        self.assertNotEqual(result.returncode, 0)
        self.assertIn("missing/stale NM-to-ENST receipt", result.stderr)
        self.assertEqual(artifact.read_bytes(), b"existing fallback artifact\n")

    def test_mapping_validation_precedes_existing_output_skip(self):
        artifact = self.folder / "motifs.pkl"
        artifact.write_bytes(b"existing artifact\n")
        self.mapping.unlink()
        result = self.run_shell('AlphaMissense_pick_missense_constrained_segments "$TEST_CONFIG"')
        self.assertNotEqual(result.returncode, 0)
        self.assertIn("Missing, empty, or unreadable nm_to_enst_tsv", result.stderr)
        self.assertEqual(artifact.read_bytes(), b"existing artifact\n")

    def test_compatible_bundle_default_is_resolved_without_config_mutation(self):
        self.config["assembly"] = "hg19"
        self.config.pop("nm_to_enst_tsv")
        before = dict(self.config)
        resolved = PREREQS.validate_mapping(self.config, ROOT)
        self.assertEqual(resolved, ROOT / PREREQS.BIOMART.RELATIVE_PATH)
        self.assertEqual(self.config, before)

    def test_mapping_install_reuses_bundle_and_registers_only_after_validation(self):
        self.config["assembly"] = "hg19"
        self.config.pop("nm_to_enst_tsv")
        result = self.run_shell('nm_to_enst_mapping_install "$TEST_CONFIG"')
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertIn("UPDATE:nm_to_enst_tsv", result.stdout)
        self.assertNotIn("UNEXPECTED_DOWNLOAD", result.stdout)

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

    def test_grch38_install_reuses_completed_bundle_without_fetching(self):
        self.config.pop("nm_to_enst_tsv")
        result = self.run_shell('''
            python() {
                if [[ $1 == */biomart_nm_to_enst.py && $2 == fetch ]]; then
                    printf 'UNEXPECTED_FETCH\\n'; return 99
                fi
                "$TEST_PYTHON" "$@"
            }
            nm_to_enst_mapping_install "$TEST_CONFIG"
        ''')
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertIn("ensembl_GRCh38_mart113_nm_to_enst.tsv", result.stderr)
        self.assertIn("UPDATE:nm_to_enst_tsv", result.stdout)
        self.assertNotIn("UNEXPECTED", result.stdout)

    def test_known_cross_assembly_defaults_are_rejected_before_acquisition(self):
        for assembly, wrong in (("hg38", "hg19"), ("hg19", "hg38")):
            with self.subTest(assembly=assembly):
                self.config["assembly"] = assembly
                self.config["nm_to_enst_tsv"] = str(ROOT / PREREQS.BIOMART.source_for_assembly(wrong)["relative_path"])
                with self.assertRaisesRegex(ValueError, "not compatible"):
                    PREREQS.nm_to_enst_path(self.config, ROOT)
                result = self.run_shell('''
                    python() {
                        if [[ $1 == */biomart_nm_to_enst.py && $2 == fetch ]]; then
                            printf 'UNEXPECTED_FETCH\\n'; return 99
                        fi
                        "$TEST_PYTHON" "$@"
                    }
                    nm_to_enst_mapping_install "$TEST_CONFIG"
                ''')
                self.assertNotEqual(result.returncode, 0)
                self.assertNotIn("UNEXPECTED", result.stdout)
                self.assertNotIn("UPDATE", result.stdout)

    def test_missing_explicit_grch38_path_does_not_fetch_the_default(self):
        self.config["nm_to_enst_tsv"] = str(self.folder / "missing_original.tsv")
        result = self.run_shell('''
            python() {
                if [[ $1 == */biomart_nm_to_enst.py && $2 == fetch ]]; then
                    printf 'UNEXPECTED_FETCH\\n'; return 99
                fi
                "$TEST_PYTHON" "$@"
            }
            nm_to_enst_mapping_install "$TEST_CONFIG"
        ''')
        self.assertNotEqual(result.returncode, 0)
        self.assertIn("Explicit nm_to_enst_tsv is missing", result.stderr)
        self.assertNotIn("UNEXPECTED", result.stdout)
        self.assertNotIn("UPDATE", result.stdout)

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

    def test_orphan_default_provenance_blocks_preflight_and_fetch_for_each_assembly(self):
        self.config.pop("nm_to_enst_tsv")
        for assembly in ("hg19", "hg38"):
            with self.subTest(assembly=assembly):
                self.config["assembly"] = assembly
                source = PREREQS.BIOMART.source_for_assembly(assembly)
                target = self.source / source["relative_path"]
                metadata = PREREQS.BIOMART.provenance_path(target)
                metadata.parent.mkdir(parents=True, exist_ok=True)
                original = '{"orphaned_original_receipt": true}\n'
                metadata.write_text(original)
                with self.assertRaisesRegex(ValueError, "existing provenance"):
                    PREREQS.validate_main(self.config, self.source, allow_missing_bundled_mapping=True)
                result = self.run_shell('''
                    BASE_DIR="$TEST_FOLDER/selected_source"
                    python() {
                        if [[ $1 == */biomart_nm_to_enst.py && $2 == fetch ]]; then
                            printf 'UNEXPECTED_FETCH\\n'; return 99
                        fi
                        "$TEST_PYTHON" "$@"
                    }
                    nm_to_enst_mapping_install "$TEST_CONFIG"
                ''')
                self.assertNotEqual(result.returncode, 0)
                self.assertIn("existing provenance", result.stderr)
                self.assertNotIn("UNEXPECTED", result.stdout)
                self.assertNotIn("UPDATE", result.stdout)
                self.assertFalse(target.exists())
                self.assertEqual(metadata.read_text(), original)

    def test_grch38_existing_table_with_stale_sidecar_is_retained_without_fetch(self):
        source = PREREQS.BIOMART.source_for_assembly("hg38")
        original = ROOT / source["relative_path"]
        target = self.source / source["relative_path"]
        target.parent.mkdir(parents=True, exist_ok=True)
        external = self.folder / "external_downloaded_grch38.tsv"
        shutil.copyfile(original, external)
        target.symlink_to(external)
        metadata = PREREQS.BIOMART.provenance_path(target)
        shutil.copyfile(PREREQS.BIOMART.provenance_path(original), metadata)
        stale = json.loads(metadata.read_text())
        stale["sha256"] = "0" * 64
        metadata.write_text(json.dumps(stale))
        before = target.read_bytes()
        self.config.pop("nm_to_enst_tsv")
        result = self.run_shell('''
            BASE_DIR="$TEST_FOLDER/selected_source"
            python() {
                if [[ $1 == */biomart_nm_to_enst.py && $2 == fetch ]]; then
                    printf 'UNEXPECTED_FETCH\\n'; return 99
                fi
                "$TEST_PYTHON" "$@"
            }
            nm_to_enst_mapping_install "$TEST_CONFIG"
        ''')
        self.assertNotEqual(result.returncode, 0)
        self.assertNotIn("UNEXPECTED", result.stdout)
        self.assertNotIn("UPDATE", result.stdout)
        self.assertEqual(target.read_bytes(), before)
        self.assertEqual(json.loads(metadata.read_text()), stale)

    def test_missing_selected_bundle_is_fetched_with_only_the_explicit_pin(self):
        self.config["assembly"] = "hg19"
        self.config["tmp_dir"] = str(self.folder / "raw-workspace")
        self.config.pop("nm_to_enst_tsv")
        self.env["TEST_COMPLETED_MAPPING"] = str(ROOT / PREREQS.BIOMART.RELATIVE_PATH)
        result = self.run_shell('''
            BASE_DIR="$TEST_FOLDER/selected_source"
            python() {
                if [[ $1 == */biomart_nm_to_enst.py && $2 == fetch ]]; then
                    "$TEST_PYTHON" -c 'import os,sys,shutil; from pathlib import Path; args=sys.argv[1:]; assert args[args.index("--source")+1]=="ensembl-grch37-mart116"; target=Path(args[args.index("--output")+1]); target.parent.mkdir(parents=True,exist_ok=True); src=os.environ["TEST_COMPLETED_MAPPING"]; shutil.copyfile(src,target); shutil.copyfile(src+".provenance.json",str(target)+".provenance.json"); print("FETCH_PIN:ensembl-grch37-mart116")' "${@:2}"
                else
                    "$TEST_PYTHON" "$@"
                fi
            }
            nm_to_enst_mapping_install "$TEST_CONFIG" "$TEST_FOLDER/raw.tsv"
        ''')
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertEqual(result.stdout.count("FETCH_PIN:ensembl-grch37-mart116"), 1)
        self.assertIn("UPDATE:nm_to_enst_tsv", result.stdout)

    def test_missing_grch38_bundle_fetches_only_release113_then_verifies(self):
        self.config["tmp_dir"] = str(self.folder / "raw-workspace")
        self.config.pop("nm_to_enst_tsv")
        source = PREREQS.BIOMART.source_for_assembly("hg38")
        self.env["TEST_COMPLETED_MAPPING"] = str(ROOT / source["relative_path"])
        result = self.run_shell('''
            BASE_DIR="$TEST_FOLDER/selected_source"
            python() {
                if [[ $1 == */biomart_nm_to_enst.py && $2 == fetch ]]; then
                    "$TEST_PYTHON" -c 'import os,sys,shutil; from pathlib import Path; args=sys.argv[1:]; assert args[args.index("--source")+1]=="ensembl-grch38-mart113"; target=Path(args[args.index("--output")+1]); assert target.name=="ensembl_GRCh38_mart113_nm_to_enst.tsv"; assert args[args.index("--raw-output")+1]==os.environ["TEST_FOLDER"]+"/raw113.tsv"; target.parent.mkdir(parents=True,exist_ok=True); src=os.environ["TEST_COMPLETED_MAPPING"]; shutil.copyfile(src,target); shutil.copyfile(src+".provenance.json",str(target)+".provenance.json"); print("FETCH_PIN:ensembl-grch38-mart113")' "${@:2}"
                else
                    "$TEST_PYTHON" "$@"
                fi
            }
            nm_to_enst_mapping_install "$TEST_CONFIG" "$TEST_FOLDER/raw113.tsv"
        ''')
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertEqual(result.stdout.count("FETCH_PIN:ensembl-grch38-mart113"), 1)
        self.assertIn("UPDATE:nm_to_enst_tsv", result.stdout)
        target = self.source / source["relative_path"]
        self.assertEqual(PREREQS.BIOMART.verify_pinned(target, "hg38")["sha256"],
                         "bb5f64c3c5eef06da92c95e50068d716c7e288629733efea6eb34dbab7bc1540")

    def test_failed_or_malformed_mapping_fetch_never_updates_configuration(self):
        self.config["assembly"] = "hg19"
        self.config["tmp_dir"] = str(self.folder / "raw-workspace")
        self.config.pop("nm_to_enst_tsv")
        result = self.run_shell('''
            BASE_DIR="$TEST_FOLDER/selected_source"
            python() {
                if [[ $1 == */biomart_nm_to_enst.py && $2 == fetch ]]; then
                    "$TEST_PYTHON" -c 'import sys; from pathlib import Path; args=sys.argv[1:]; target=Path(args[args.index("--output")+1]); target.parent.mkdir(parents=True,exist_ok=True); target.write_text("failed HTTP response\\n")' "${@:2}"
                else
                    "$TEST_PYTHON" "$@"
                fi
            }
            nm_to_enst_mapping_install "$TEST_CONFIG" "$TEST_FOLDER/raw.tsv"
        ''')
        self.assertNotEqual(result.returncode, 0)
        self.assertNotIn("UPDATE", result.stdout)
        self.assertEqual(json.loads(self.config_file.read_text()), self.config)

    def test_matching_mapping_anchor_and_output_receipts_allow_existing_cache_skip(self):
        for function, name in (("AlphaMissense_pick_missense_constrained_segments", "motifs.pkl"),
                               ("prepare_pm1_regions", "pm1_regions.hg38.pkl")):
            with self.subTest(function=function):
                artifact = self.folder / name
                artifact.write_bytes(b"fixture completed scientific output\n")
                os.utime(artifact, (artifact.stat().st_atime + 10, artifact.stat().st_mtime + 10))
                anchors = [ROOT / "data/gnomAD_RMC/gnomad_v2.1.1_rmc_flat.hg38.tsv",
                           ROOT / "data/HCSeeker/HC_spots.hg38.tsv"]
                identity = PREREQS.BIOMART.mapping_identity(self.mapping, anchors)
                PREREQS.BIOMART.atomic_json(Path(str(artifact) + ".nm_to_enst.json"),
                    {"mapping_identity": identity, "artifact_bytes": artifact.stat().st_size,
                     "artifact_sha256": PREREQS.BIOMART.sha256(artifact)})
                result = self.run_shell('''
                    python() {
                        if [[ $1 == */installer_prerequisites.py || $1 == */biomart_nm_to_enst.py ]]; then
                            "$TEST_PYTHON" "$@"
                        else
                            printf 'UNEXPECTED_SCIENTIFIC_COMMAND\\n'; return 99
                        fi
                    }
                ''' + f'{function} "$TEST_CONFIG"')
                self.assertEqual(result.returncode, 0, result.stderr)
                self.assertIn("skip this step", result.stderr)
                self.assertNotIn("UNEXPECTED", result.stdout)

    def test_gofcards_hg38_install_passes_distinct_original_references(self):
        # Select an external final output so this exercises the normalized/build route.
        self.config["gofcards_exact_gof_cache"] = str(self.folder / "external_final_gof.json.gz")
        result = self.run_shell('''
            vep() { printf 'ensembl-vep 113.0\\n'; }
            bcftools() { :; }
            CrossMap() { :; }
            validate_gofcards_exact_gof_cache() { return 0; }
            python() {
                if [[ $1 == */build_gofcards_exact_gof_cache.py ]]; then
                    "$TEST_PYTHON" -c 'import json,sys; print("BUILD:"+json.dumps(sys.argv[1:]))' "$@"
                else
                    "$TEST_PYTHON" "$@"
                fi
            }
            gofcards_exact_gof_cache_install "$TEST_CONFIG"
        ''')
        self.assertEqual(result.returncode, 0, result.stderr)
        call = next(json.loads(line[6:]) for line in result.stdout.splitlines() if line.startswith("BUILD:"))
        self.assertEqual(call[call.index("--hg19-fasta") + 1], str(self.references["hg19"]))
        self.assertEqual(call[call.index("--hg38-fasta") + 1], self.config["ref_genome"])
        self.assertEqual(call[call.index("--vep-cache-version") + 1], "113")
        self.assertNotIn("UNEXPECTED_DOWNLOAD", result.stdout)

    def test_gofcards_missing_opposite_reference_stops_before_work(self):
        # Select an external final output so this exercises the normalized/build route.
        self.config["gofcards_exact_gof_cache"] = str(self.folder / "external_final_gof.json.gz")
        self.references["hg19"].unlink()
        result = self.run_shell('''
            priva_verify_vep_release() { [[ $1 == 113 ]]; }
            gofcards_exact_gof_cache_install "$TEST_CONFIG" --main-preflight
        ''')
        self.assertNotEqual(result.returncode, 0)
        self.assertIn("gofcards_hg19_fasta", result.stderr)
        self.assertFalse((self.folder / "build_work").exists())
        self.assertNotIn("UNEXPECTED_DOWNLOAD", result.stdout)

    def test_current_normalized_cache_reuses_without_unused_build_inputs(self):
        # Select an external final output so this exercises the normalized/build route.
        self.config["gofcards_exact_gof_cache"] = str(self.folder / "external_final_gof.json.gz")
        work = self.folder / "build_work"
        work.mkdir()
        target = work / "gofcards_exact_gof.normalized.json.gz"
        target.write_text("validated by the shell fixture validator\n")
        self.config["gofcards_builder_script"] = str(self.folder / "absent_builder.py")
        self.config["gofcards_mechanism_review_tsv"] = str(self.folder / "absent_review.tsv")
        self.config.pop("gofcards_hg19_fasta")
        self.config.pop("gofcards_hg38_fasta")
        # Ensure no unchanged source predates the normalized cache accidentally.
        os.utime(target, (target.stat().st_atime + 5, target.stat().st_mtime + 5))
        result = self.run_shell('''
            validate_gofcards_exact_gof_cache() { return 0; }
            gofcards_exact_gof_cache_install "$TEST_CONFIG"
        ''')
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertIn("normalized cache is current", result.stderr)
        self.assertNotIn("UNEXPECTED_DOWNLOAD", result.stdout)


if __name__ == "__main__":
    unittest.main()
