"""Installer dependency/command contracts without downloads or scientific runs."""
import os
import json
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest

ROOT = Path(__file__).resolve().parents[1]
INSTALLER = ROOT / "scripts/install_utils.sh"


class InstallerResources(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory(prefix="priva_resources_", dir=os.environ.get("PRIVA_TEST_TMPDIR", ROOT / "tests"))
        self.addCleanup(self.temp.cleanup)
        self.folder = Path(self.temp.name)
        self.reference = self.folder / "reference.fa"
        self.reference.write_text(">test\nA\n")
        self.chain = self.folder / "test.chain.gz"
        self.chain.write_text("stub chain\n")
        self.env = dict(os.environ, TEST_FOLDER=str(self.folder), TEST_REFERENCE=str(self.reference),
                        TEST_CHAIN=str(self.chain))

    def run_shell(self, body):
        return subprocess.run(["bash", "-c", 'source "$1"; log() { printf "%s\\n" "$*" >&2; }; '
                               'wget() { printf "UNEXPECTED_NETWORK\\n"; return 98; }; ' + body,
                               "test", str(INSTALLER)], text=True, capture_output=True, env=self.env)

    def gnomad_settings(self):
        return '''
            read_yaml() {
                case "$2" in
                    assembly) printf '%s\\n' "${TEST_ASSEMBLY:-GRCh37}" ;;
                    ref_genome) printf '%s\\n' "$TEST_REFERENCE" ;;
                    gnomad_vcf_dir) printf '%s\\n' "$TEST_FOLDER" ;;
                esac
            }
            update_yaml() { printf 'UPDATE:%s\\n' "$*"; }
        '''

    def test_liftover_receives_configured_reference_without_assembly_in_filename(self):
        result = self.run_shell(self.gnomad_settings() + '''
            check_vcf_validity() { return 0; }
            gnomAD_liftover() { printf 'LIFT:%s\\n' "$3"; }
            parallel() { return 0; }
            gnomAD_install config "$TEST_REFERENCE"
        ''')
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertIn("LIFT:" + str(self.reference), result.stdout)
        self.assertIn("gnomad.joint.v4.1.sites.hg19.chrX.vcf.gz", result.stdout)

    def test_missing_grch37_reference_stops_before_resource_access(self):
        self.env["TEST_REFERENCE"] = str(self.folder / "absent.fa")
        result = self.run_shell(self.gnomad_settings() + '''
            check_vcf_validity() { printf 'RESOURCE_ACCESS\\n'; return 0; }
            gnomAD_install config
        ''')
        self.assertNotEqual(result.returncode, 0)
        self.assertNotIn("RESOURCE_ACCESS", result.stdout)

    def test_failed_gnomad_download_stops_without_updating_config(self):
        self.env["TEST_ASSEMBLY"] = "GRCh38"
        result = self.run_shell(self.gnomad_settings() + '''
            check_vcf_validity() { return 1; }
            wget() { printf 'DOWNLOAD\\n'; return 71; }
            gnomAD_install config
        ''')
        self.assertNotEqual(result.returncode, 0)
        self.assertEqual(result.stdout.count("DOWNLOAD"), 1)
        self.assertNotIn("UPDATE", result.stdout)

    def test_all_chromosomes_are_scheduled_and_failure_stops_concatenation(self):
        result = self.run_shell('''
            mkdir() { return 0; }
            parallel() {
                printf 'PARALLEL:%s\\n' "$*"
                [[ "$*" == *--dry-run* ]] && return 0
                return 72
            }
            bcftools_concatvcfs() { printf 'CONCAT_CALLED\\n'; }
            gnomAD_liftover "$TEST_FOLDER/gnomad.chr1.vcf.bgz" "$TEST_CHAIN" "$TEST_REFERENCE" 4
        ''')
        self.assertNotEqual(result.returncode, 0)
        calls = [line.split(" ::: ")[1].split() for line in result.stdout.splitlines() if line.startswith("PARALLEL:")]
        expected = ["chr" + str(i) for i in range(1, 23)] + ["chrX", "chrY"]
        self.assertEqual(calls, [expected, expected])
        self.assertNotIn("CONCAT_CALLED", result.stdout)
        self.assertIn("refusing partial concatenation", result.stderr)

    def test_invalid_gnomad_download_never_updates_configuration(self):
        self.env["TEST_ASSEMBLY"] = "GRCh38"
        result = self.run_shell(self.gnomad_settings() + '''
            check_vcf_validity() { return 1; }
            wget() { printf 'DOWNLOAD\\n'; }
            gnomAD_install config
        ''')
        self.assertNotEqual(result.returncode, 0)
        self.assertEqual(result.stdout.count("DOWNLOAD"), 2)
        self.assertNotIn("UPDATE", result.stdout)
        self.assertIn("is invalid", result.stderr)

    def test_failed_gnomad_reindex_never_updates_configuration(self):
        source = self.folder / "gnomad.joint.v4.1.sites.chr1.vcf.bgz"
        source.write_text("stub VCF\n")
        index = Path(str(source) + ".tbi")
        index.write_text("stub index\n")
        os.utime(index, (1, 1))
        self.env["TEST_ASSEMBLY"] = "GRCh38"
        result = self.run_shell(self.gnomad_settings() + '''
            check_vcf_validity() { return 0; }
            tabix() { printf 'INDEX_FAILED\\n'; return 76; }
            gnomAD_install config
        ''')
        self.assertNotEqual(result.returncode, 0)
        self.assertIn("INDEX_FAILED", result.stdout)
        self.assertNotIn("UPDATE", result.stdout)
        self.assertNotIn("UNEXPECTED_NETWORK", result.stdout)

    def test_missing_liftover_chain_stops_before_scheduling(self):
        self.env["TEST_CHAIN"] = str(self.folder / "missing.chain.gz")
        result = self.run_shell('''
            parallel() { printf 'UNEXPECTED_PARALLEL\\n'; }
            gnomAD_liftover "$TEST_FOLDER/gnomad.chr1.vcf.bgz" "$TEST_CHAIN" "$TEST_REFERENCE" 4
        ''')
        self.assertNotEqual(result.returncode, 0)
        self.assertNotIn("UNEXPECTED_PARALLEL", result.stdout)

    def test_failed_lifted_vcf_query_stops_before_splitting(self):
        result = self.run_shell('''
            check_vcf_validity() { [[ "$1" != *.hg19.vcf.gz ]]; }
            crossmap_liftover_hg382hg19() { return 0; }
            bcftools() { printf 'chr1\\n'; return 77; }
            display_vcf() { printf 'UNEXPECTED_DISPLAY\\n'; }
            gnomAD_liftover_per_chromosome "$TEST_FOLDER/gnomad.chr1.vcf.bgz" "$TEST_FOLDER" "$TEST_REFERENCE" "$TEST_CHAIN"
        ''')
        self.assertNotEqual(result.returncode, 0)
        self.assertIn("Failed to read chromosomes", result.stderr)
        self.assertNotIn("UNEXPECTED_DISPLAY", result.stdout)

    def test_failed_lifted_vcf_split_or_index_propagates(self):
        for stage in ("view", "index"):
            with self.subTest(stage=stage):
                self.env["TEST_FAIL_STAGE"] = stage
                result = self.run_shell('''
                    check_vcf_validity() { [[ "$1" != *.hg19.vcf.gz ]]; }
                    crossmap_liftover_hg382hg19() { return 0; }
                    bcftools() {
                        if [[ "$1" == query ]]; then printf 'chr1\\n'; return 0; fi
                        [[ "$1" != "$TEST_FAIL_STAGE" ]]
                    }
                    display_vcf() { printf 'UNEXPECTED_DISPLAY\\n'; }
                    gnomAD_liftover_per_chromosome "$TEST_FOLDER/gnomad.chr1.vcf.bgz" "$TEST_FOLDER" "$TEST_REFERENCE" "$TEST_CHAIN"
                ''')
                self.assertNotEqual(result.returncode, 0)
                self.assertIn("Failed to split or index", result.stderr)
                self.assertNotIn("UNEXPECTED_DISPLAY", result.stdout)

    def test_unsupported_gnomad_assembly_stops_before_resource_access(self):
        self.env["TEST_ASSEMBLY"] = "GRCh39"
        result = self.run_shell(self.gnomad_settings() + '''
            check_vcf_validity() { printf 'UNEXPECTED_RESOURCE\\n'; }
            gnomAD_install config
        ''')
        self.assertNotEqual(result.returncode, 0)
        self.assertNotIn("UNEXPECTED", result.stdout)

    def plugin_stubs(self):
        return '''
            CONDA_PREFIX="$TEST_FOLDER/envs/priva_acmg"
            for name in UTRAnnotator_install LOEUF_install AlphaMissense_install AlphaMissense_pick_missense_constrained_segments SpliceAI_install SpliceVault_install MaveDB_install PrimateAI_install; do
                eval "$name() { :; }"
            done
            AlphaMissense_anno() { printf 'ANNOTATION\\n'; }
        '''

    def test_conservation_is_ready_before_alphamissense_annotation(self):
        result = self.run_shell(self.plugin_stubs() + '''
            Conservation_install() { printf 'CONSERVATION\\n'; }
            VEP_plugins_install plugins caches hg19 config priva_acmg
        ''')
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertEqual(result.stdout.splitlines(), ["CONSERVATION", "ANNOTATION"])

    def test_conservation_failure_prevents_alphamissense_annotation(self):
        result = self.run_shell(self.plugin_stubs() + '''
            Conservation_install() { return 73; }
            VEP_plugins_install plugins caches hg19 config priva_acmg
        ''')
        self.assertNotEqual(result.returncode, 0)
        self.assertNotIn("ANNOTATION", result.stdout)

    def test_existing_conservation_file_routing_is_reused_for_both_assemblies(self):
        # Routing fixture only; dedicated native-BigWig tests exercise validation.
        for assembly in ("hg19", "GRCh38"):
            with self.subTest(assembly=assembly):
                self.env["TEST_ASSEMBLY"] = assembly
                result = self.run_shell('''
                    read_yaml() { printf '%s\\n' "$TEST_REFERENCE"; }
                    python() { printf 'VALIDATE:%s\\n' "$*"; }
                    wget() { printf 'UNEXPECTED_DOWNLOAD\\n'; return 74; }
                    update_yaml() { printf 'UNEXPECTED_UPDATE\\n'; }
                    Conservation_install config "$TEST_FOLDER" "$TEST_ASSEMBLY"
                ''')
                self.assertEqual(result.returncode, 0, result.stderr)
                self.assertNotIn("UNEXPECTED", result.stdout)
                self.assertIn("conservation_resource_check.py --input " + str(self.reference), result.stdout)

    def test_conservation_download_failure_never_updates_configuration(self):
        result = self.run_shell('''
            read_yaml() { case "$2" in conservation_file) printf '%s/Conservation/gerp_conservation_scores.homo_sapiens.GRCh37.bw\\n' "$TEST_FOLDER" ;; ref_genome) printf '%s\\n' "$TEST_REFERENCE" ;; esac; }
            wget() { return 75; }
            update_yaml() { printf 'UNEXPECTED_UPDATE\\n'; }
            Conservation_install config "$TEST_FOLDER" hg19
        ''')
        self.assertNotEqual(result.returncode, 0)
        self.assertNotIn("UNEXPECTED_UPDATE", result.stdout)
        self.assertFalse((self.folder / "Conservation/gerp_conservation_scores.homo_sapiens.GRCh37.bw").exists())

    def test_grch38_conservation_default_keeps_archived91_mammals_release114(self):
        result = self.run_shell('''
            read_yaml() { case "$2" in conservation_file) printf '%s/Conservation/gerp_conservation_scores.homo_sapiens.GRCh38.bw\\n' "$TEST_FOLDER" ;; conservation_GRCh38_url) printf 'null\\n' ;; ref_genome) printf '%s\\n' "$TEST_REFERENCE" ;; esac; }
            wget() { printf 'DOWNLOAD:%s\\n' "$2"; return 75; }
            update_yaml() { printf 'UNEXPECTED_UPDATE\\n'; }
            Conservation_install config "$TEST_FOLDER" hg38
        ''')
        self.assertNotEqual(result.returncode, 0)
        self.assertIn("DOWNLOAD:https://ftp.ensembl.org/pub/release-114/compara/conservation_scores/91_mammals.gerp_conservation_score/gerp_conservation_scores.homo_sapiens.GRCh38.bw", result.stdout)
        self.assertNotIn("current_compara", result.stdout)
        self.assertNotIn("UNEXPECTED_UPDATE", result.stdout)

    def test_grch38_conservation_uses_explicit_configured_source_without_substitution(self):
        result = self.run_shell('''
            read_yaml() { case "$2" in conservation_file) printf '%s/Conservation/gerp_conservation_scores.homo_sapiens.GRCh38.bw\\n' "$TEST_FOLDER" ;; conservation_GRCh38_url) printf 'https://example.invalid/explicit-pinned91.bw\\n' ;; ref_genome) printf '%s\\n' "$TEST_REFERENCE" ;; esac; }
            wget() { printf 'DOWNLOAD:%s\\n' "$2"; return 75; }
            Conservation_install config "$TEST_FOLDER" GRCh38
        ''')
        self.assertNotEqual(result.returncode, 0)
        self.assertIn("DOWNLOAD:https://example.invalid/explicit-pinned91.bw", result.stdout)
        self.assertNotIn("release-114", result.stdout)

    def test_empty_conservation_download_never_updates_configuration(self):
        result = self.run_shell('''
            read_yaml() { case "$2" in conservation_file) printf '%s/Conservation/gerp_conservation_scores.homo_sapiens.GRCh37.bw\\n' "$TEST_FOLDER" ;; ref_genome) printf '%s\\n' "$TEST_REFERENCE" ;; esac; }
            wget() { return 0; }
            update_yaml() { printf 'UNEXPECTED_UPDATE\\n'; }
            Conservation_install config "$TEST_FOLDER" hg19
        ''')
        self.assertNotEqual(result.returncode, 0)
        self.assertNotIn("UNEXPECTED_UPDATE", result.stdout)
        self.assertFalse((self.folder / "Conservation/gerp_conservation_scores.homo_sapiens.GRCh37.bw").exists())

    def test_loftee_complete_resources_reuse_under_bash_nounset(self):
        # Routing only: native gzip/FAI/SQLite/BigWig validation has dedicated tests.
        Path(str(self.reference) + ".fai").write_text("fixture index\n")
        plugins = self.folder / "plugins"
        caches = self.folder / "plugin_cache"
        caches.mkdir()
        for assembly in ("hg19", "GRCh37", "hg38", "GRCh38"):
            with self.subTest(assembly=assembly):
                repository = plugins / "LoFtee_repos" / f"loftee-{assembly}"
                repository.mkdir(parents=True)
                (repository / "LoF.pm").write_text("# fixture plugin\n")
                body = r'''
                    source "$1"
                    log() { printf '%s\n' "$*" >&2; }
                    read_yaml() {
                        case "$2" in
                            assembly) printf '%s\n' "$TEST_ASSEMBLY" ;;
                            vep_plugins_dir) printf '%s/plugins\n' "$TEST_FOLDER" ;;
                            human_ancestor_fasta|loftee_conservation_file|gerp_bigwig) printf '%s\n' "$TEST_REFERENCE" ;;
                            *) return 96 ;;
                        esac
                    }
                    git() { printf 'UNEXPECTED_GIT\n'; return 97; }
                    wget() { printf 'UNEXPECTED_DOWNLOAD\n'; return 98; }
                    update_yaml() { printf 'UNEXPECTED_UPDATE\n'; return 99; }
                    validate_loftee_resources() { return 0; }
                    unset loftee_cache_dir
                    LoFtee_install "$TEST_FOLDER/plugin_cache" config
                '''
                environment = dict(self.env, TEST_ASSEMBLY=assembly, TMPDIR=str(self.folder))
                result = subprocess.run(["bash", "--noprofile", "--norc", "-u", "-c", body,
                                         "test", str(INSTALLER)], env=environment, text=True,
                                        capture_output=True, timeout=30)
                self.assertEqual(result.returncode, 0, result.stderr)
                self.assertIn("already installed", result.stderr)
                self.assertNotIn("UNEXPECTED", result.stdout)
                self.assertNotIn("missing operand", result.stderr)
                self.assertNotIn("unbound variable", result.stderr)
                self.assertTrue((caches / "LoFtee" / assembly).is_dir())

    def test_alphamissense_motif_output_stays_at_configured_cluster_path(self):
        (self.folder / "original_mapping.tsv").write_text(
            "Gene stable ID\tTranscript stable ID\tHGNC symbol\tUniProtKB/Swiss-Prot ID\tRefSeq mRNA ID\n"
            "ENSG000001\tENST000001\tGENE\tP00001\tNM_000001\n"
        )
        config = self.folder / "mapping_config.json"
        config.write_text(json.dumps({"assembly": "hg19", "nm_to_enst_tsv": str(self.folder / "original_mapping.tsv")}))
        self.env["TEST_MAPPING_CONFIG"] = str(config)
        self.env["TEST_PYTHON"] = sys.executable
        result = self.run_shell('''
            read_yaml() {
                case "$2" in
                    assembly) printf 'hg19\\n' ;;
                    alphamissense_intolerant_motifs) printf '%s/motifs.pkl\\n' "$TEST_FOLDER" ;;
                    alphamissense_vep_vcf) printf '%s/input.vcf.gz\\n' "$TEST_FOLDER" ;;
                    nm_to_enst_tsv) printf '%s/original_mapping.tsv\\n' "$TEST_FOLDER" ;;
                    threads) printf '2\\n' ;;
                esac
            }
            python() {
                case "$1" in
                    */installer_prerequisites.py|*/biomart_nm_to_enst.py) "$TEST_PYTHON" "$@" ;;
                    */alphamissense_annotation_base.py)
                        case "$2" in
                            path) printf '%s/input.vcf.gz\\n' "$TEST_FOLDER" ;;
                            check) return 0 ;;
                            motif-check) return 1 ;; # Exercise fresh generation.
                            *) return 97 ;;
                        esac ;;
                    */mapped_artifact_publish.py)
                        [[ $2 != recover ]] || return 0
                        [[ $2 == run ]] || return 98
                        printf 'PYTHON:%s\\n' "$*"
                        "$TEST_PYTHON" -c 'import sys; from pathlib import Path; args=sys.argv[1:]; Path(args[args.index("--output")+1]).write_text("fixture output\\n")' "$@" ;;
                    *) return 99 ;;
                esac
            }
            update_yaml() { printf 'UPDATE:%s\\n' "$*"; }
            AlphaMissense_pick_missense_constrained_segments "$TEST_MAPPING_CONFIG"
        ''')
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertIn("--output " + str(self.folder / "motifs.pkl"), result.stdout)
        self.assertNotIn("/paedyl01", result.stdout + result.stderr)


if __name__ == "__main__":
    unittest.main()
