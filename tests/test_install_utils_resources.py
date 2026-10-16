"""Installer dependency/command contracts without downloads or scientific runs."""
import os
from pathlib import Path
import subprocess
import tempfile
import unittest

ROOT = Path(__file__).resolve().parents[1]
INSTALLER = ROOT / "scripts/install_utils.sh"


class InstallerResources(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory(prefix="priva_resources_", dir=ROOT / "tests")
        self.addCleanup(self.temp.cleanup)
        self.folder = Path(self.temp.name)
        self.reference = self.folder / "reference.fa"
        self.reference.write_text(">test\nA\n")
        self.env = dict(os.environ, TEST_FOLDER=str(self.folder), TEST_REFERENCE=str(self.reference))

    def run_shell(self, body):
        return subprocess.run(["bash", "-c", 'source "$1"; log() { printf "%s\\n" "$*" >&2; }; ' + body,
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
            gnomAD_liftover "$TEST_FOLDER/gnomad.chr1.vcf.bgz" chain "$TEST_REFERENCE" 4
        ''')
        self.assertNotEqual(result.returncode, 0)
        calls = [line.split(" ::: ")[1].split() for line in result.stdout.splitlines() if line.startswith("PARALLEL:")]
        expected = ["chr" + str(i) for i in range(1, 23)] + ["chrX", "chrY"]
        self.assertEqual(calls, [expected, expected])
        self.assertNotIn("CONCAT_CALLED", result.stdout)
        self.assertIn("refusing partial concatenation", result.stderr)

    def plugin_stubs(self):
        return '''
            CONDA_PREFIX="$TEST_FOLDER/envs/priva_acmg"
            for name in UTRAnnotator_install LOEUF_install AlphaMissense_install AlphaMissense_pick_missense_constrained_segments SpliceAI_profile_install SpliceVault_install MaveDB_install PrimateAI_install; do
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

    def test_existing_conservation_file_is_reused_for_both_assemblies(self):
        for assembly in ("hg19", "GRCh38"):
            with self.subTest(assembly=assembly):
                self.env["TEST_ASSEMBLY"] = assembly
                result = self.run_shell('''
                    read_yaml() { printf '%s\\n' "$TEST_REFERENCE"; }
                    wget() { printf 'UNEXPECTED_DOWNLOAD\\n'; return 74; }
                    update_yaml() { printf 'UNEXPECTED_UPDATE\\n'; }
                    Conservation_install config "$TEST_FOLDER" "$TEST_ASSEMBLY"
                ''')
                self.assertEqual(result.returncode, 0, result.stderr)
                self.assertNotIn("UNEXPECTED", result.stdout)

    def test_conservation_download_failure_never_updates_configuration(self):
        result = self.run_shell('''
            read_yaml() { printf '%s\\n' "$TEST_FOLDER/absent.bw"; }
            wget() { return 75; }
            update_yaml() { printf 'UNEXPECTED_UPDATE\\n'; }
            Conservation_install config "$TEST_FOLDER" hg19
        ''')
        self.assertNotEqual(result.returncode, 0)
        self.assertNotIn("UNEXPECTED_UPDATE", result.stdout)
        self.assertFalse((self.folder / "Conservation/gerp_conservation_scores.homo_sapiens.GRCh37.bw").exists())

    def test_alphamissense_motif_output_stays_at_configured_cluster_path(self):
        result = self.run_shell('''
            read_yaml() {
                case "$2" in
                    assembly) printf 'hg19\\n' ;;
                    alphamissense_intolerant_motifs) printf '%s/motifs.pkl\\n' "$TEST_FOLDER" ;;
                    alphamissense_vep_vcf) printf '%s/input.vcf.gz\\n' "$TEST_FOLDER" ;;
                    threads) printf '2\\n' ;;
                esac
            }
            python() { printf 'PYTHON:%s\\n' "$*"; }
            update_yaml() { printf 'UPDATE:%s\\n' "$*"; }
            AlphaMissense_pick_missense_constrained_segments config
        ''')
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertIn("--output " + str(self.folder / "motifs.pkl"), result.stdout)
        self.assertNotIn("/paedyl01", result.stdout + result.stderr)


if __name__ == "__main__":
    unittest.main()
