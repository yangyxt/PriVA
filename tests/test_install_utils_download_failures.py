"""Failure propagation for original installer downloads; no network or scientific runs."""
import gzip
import hashlib
import os
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest

ROOT = Path(__file__).resolve().parents[1]
INSTALLER = ROOT / "scripts/install_utils.sh"


class InstallerDownloadFailures(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory(prefix="priva_download_failures_",
                                               dir=os.environ.get("PRIVA_TEST_TMPDIR", os.environ.get("TMPDIR")))
        self.addCleanup(self.temp.cleanup)
        self.folder = Path(self.temp.name)
        self.cache = self.folder / "cache"
        self.cache.mkdir()
        (self.cache / "InterPro").mkdir()
        (self.cache / "ClinGen").mkdir()
        self.data = self.folder / "data"
        (self.data / "InterPro").mkdir(parents=True)
        (self.data / "clingen").mkdir()
        self.xml_fixture = self.folder / "fixture.xml.gz"
        with gzip.open(self.xml_fixture, "wb") as stream:
            stream.write(b"<interprodb/>\n")
        self.env = dict(os.environ, TEST_FOLDER=str(self.folder), TEST_CACHE=str(self.cache),
                        TEST_DATA=str(self.data), TEST_ASSEMBLY="hg19", TEST_PYTHON=sys.executable,
                        TEST_INTERPRO_XML=str(self.xml_fixture),
                        TEST_INTERPRO_MAP=str(self.cache / "InterPro/Interpro_entry_mapping.pkl.gz"),
                        TEST_UTR_FILE="")

    def run_shell(self, body):
        settings = r'''
            DATA_DIR="$TEST_DATA"
            log() { printf '%s\n' "$*" >&2; }
            read_yaml() {
                case "$2" in
                    assembly) printf '%s\n' "$TEST_ASSEMBLY" ;;
                    utr_annotator_file) printf '%s\n' "$TEST_UTR_FILE" ;;
                    interpro_mapping_pickle) printf '%s\n' "$TEST_INTERPRO_MAP" ;;
                    clingen_map) printf '%s/ClinGen/clingen_map.hg19.pkl.gz\n' "$TEST_CACHE" ;;
                    gene_dosage_sensitivity) printf '%s/ClinGen/gene_dosage_sensitivity.hg19.tsv\n' "$TEST_CACHE" ;;
                    interpro_xml_gz_url) printf 'https://fixtures.invalid/interpro.xml.gz\n' ;;
                    alphamissense_vep_vcf) printf '%s/alpha.vcf.gz\n' "$TEST_FOLDER" ;;
                    clinvar_vcf) printf '%s/clinvar.vcf.gz\n' "$TEST_FOLDER" ;;
                esac
            }
            update_yaml() { printf 'UPDATE:%s %s\n' "$2" "$3"; }
            wget() { printf 'UNEXPECTED_NETWORK\n'; return 98; }
            git() { printf 'UNEXPECTED_GIT\n'; return 97; }
            python() { printf 'UNEXPECTED_PYTHON\n'; return 96; }
        '''
        return subprocess.run(["bash", "-c", 'source "$1"; ' + settings + body,
                               "test", str(INSTALLER)], text=True, capture_output=True, env=self.env)

    def clinvar_inputs(self, dosage=True):
        (self.cache / "ClinGen/clingen_erepo.tsv").write_text("stub eRepo\n")
        if dosage:
            (self.cache / "ClinGen/gene_dosage_sensitivity.hg19.tsv").write_text(
                "# comment\n#Gene\tScore\nGENE\t3\n")

    def test_failed_utr_clone_never_updates_configuration(self):
        result = self.run_shell(r'''
            git() { printf 'CLONE_FAILED\n'; return 71; }
            UTRAnnotator_install config "$TEST_CACHE"
        ''')
        self.assertNotEqual(result.returncode, 0)
        self.assertIn("CLONE_FAILED", result.stdout)
        self.assertNotIn("UPDATE:", result.stdout)

    def test_missing_utr_resource_after_clone_never_updates_configuration(self):
        result = self.run_shell(r'''
            git() { return 0; }
            UTRAnnotator_install config "$TEST_CACHE"
        ''')
        self.assertNotEqual(result.returncode, 0)
        self.assertIn("resource is missing or empty", result.stderr)
        self.assertNotIn("UPDATE:", result.stdout)

    def test_existing_nonempty_utr_resource_is_reused(self):
        resource = self.folder / "configured_utr.txt"
        resource.write_text("stub UTR resource\n")
        self.env["TEST_UTR_FILE"] = str(resource)
        result = self.run_shell('UTRAnnotator_install config "$TEST_CACHE"')
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertNotIn("UNEXPECTED", result.stdout)
        self.assertNotIn("UPDATE:", result.stdout)

    def test_utr_clone_keeps_original_source_and_assembly_selection(self):
        for assembly, resource_assembly in (("hg19", "GRCh37"), ("hg38", "GRCh38")):
            with self.subTest(assembly=assembly):
                self.env["TEST_ASSEMBLY"] = assembly
                result = self.run_shell(r'''
                    git() {
                        printf 'GIT:%s\n' "$*"
                        mkdir -p UTRannotator
                        printf 'stub hg19\n' > UTRannotator/uORF_5UTR_GRCh37_PUBLIC.txt
                        printf 'stub hg38\n' > UTRannotator/uORF_5UTR_GRCh38_PUBLIC.txt
                    }
                    UTRAnnotator_install config "$TEST_CACHE"
                ''')
                self.assertEqual(result.returncode, 0, result.stderr)
                target = self.cache / ("UTRannotator/uORF_5UTR_" + resource_assembly + "_PUBLIC.txt")
                self.assertIn("UPDATE:utr_annotator_file " + str(target), result.stdout)
                if assembly == "hg19":
                    self.assertIn("GIT:clone https://github.com/Ensembl/UTRannotator", result.stdout)

    def test_failed_interpro_download_stops_before_parser(self):
        result = self.run_shell(r'''
            wget() { printf 'DOWNLOAD_FAILED\n'; return 72; }
            InterPro_parsing config
        ''')
        self.assertNotEqual(result.returncode, 0)
        self.assertIn("DOWNLOAD_FAILED", result.stdout)
        self.assertNotIn("UNEXPECTED_PYTHON", result.stdout)
        self.assertNotIn("UPDATE:", result.stdout)

    def test_interpro_failed_transfer_resumes_without_publishing_incomplete_source(self):
        partial = self.cache / "InterPro/interpro.xml.gz.partial.gz"
        target = self.cache / "InterPro/interpro.xml.gz"
        failed = self.run_shell(r'''
            wget() { printf 'unfinished source bytes' > "${@: -1}"; return 72; }
            InterPro_parsing config
        ''')
        self.assertNotEqual(failed.returncode, 0)
        self.assertEqual(partial.read_bytes(), b"unfinished source bytes")
        self.assertFalse(target.exists())
        resumed = self.run_shell(r'''
            wget() {
                [[ "$1" == -c && -s "${@: -1}" ]] || return 75
                cp "$TEST_INTERPRO_XML" "${@: -1}"
            }
            python() { printf 'PARSER_CALLED\n'; return 73; }
            InterPro_parsing config
        ''')
        self.assertIn("PARSER_CALLED", resumed.stdout)
        self.assertEqual(target.read_bytes(), self.xml_fixture.read_bytes())
        self.assertFalse(partial.exists())
        self.assertNotIn("UPDATE:", resumed.stdout)

    def test_invalid_interpro_download_is_never_published(self):
        result = self.run_shell(r'''
            wget() { printf 'invalid gzip' > "${@: -1}"; }
            InterPro_parsing config
        ''')
        self.assertNotEqual(result.returncode, 0)
        self.assertIn("Invalid downloaded InterPro gzip", result.stderr)
        self.assertFalse((self.cache / "InterPro/interpro.xml.gz").exists())
        self.assertEqual((self.cache / "InterPro/interpro.xml.gz.partial.gz").read_bytes(), b"invalid gzip")
        self.assertNotIn("UNEXPECTED_PYTHON", result.stdout)
        self.assertNotIn("UPDATE:", result.stdout)

    def test_empty_interpro_download_stops_before_parser(self):
        result = self.run_shell(r'''
            wget() { return 0; }
            InterPro_parsing config
        ''')
        self.assertNotEqual(result.returncode, 0)
        self.assertIn("resource is missing or empty", result.stderr)
        self.assertNotIn("UNEXPECTED_PYTHON", result.stdout)
        self.assertNotIn("UPDATE:", result.stdout)

    def test_interpro_parser_failure_never_updates_configuration(self):
        result = self.run_shell(r'''
            wget() { cp "$TEST_INTERPRO_XML" "${@: -1}"; }
            python() { printf 'PARSER_FAILED\n'; return 73; }
            InterPro_parsing config
        ''')
        self.assertNotEqual(result.returncode, 0)
        self.assertIn("PARSER_FAILED", result.stdout)
        self.assertNotIn("UPDATE:", result.stdout)
        self.assertTrue((self.cache / "InterPro/interpro.xml.gz").is_file())

    def test_interpro_success_keeps_original_parser_arguments(self):
        result = self.run_shell(r'''
            wget() { printf 'DOWNLOAD:%s\n' "$*"; cp "$TEST_INTERPRO_XML" "${@: -1}"; }
            python() {
                printf 'PARSER:%s\n' "$*"
                while [[ $# -gt 0 ]]; do
                    if [[ "$1" == --mapping_output ]]; then printf 'stub mapping\n' > "$2"; break; fi
                    shift
                done
            }
            InterPro_parsing config
        ''')
        self.assertEqual(result.returncode, 0, result.stderr)
        interpro = self.cache / "InterPro"
        report = interpro / "Interpro_format_report.txt"
        self.assertIn("DOWNLOAD:-c https://fixtures.invalid/interpro.xml.gz -O " +
                      str(interpro / "interpro.xml.gz.partial.gz"), result.stdout)
        self.assertIn("protein_domain_mapping.py --xml " + str(interpro / "interpro.xml.gz") +
                      " --mapping_output " + str(interpro / "Interpro_entry_mapping.pkl.gz") +
                      " --vcf " + str(self.folder / "alpha.vcf.gz") +
                      " --report " + str(report), result.stdout)
        self.assertIn("UPDATE:interpro_mapping_pickle " +
                      str(interpro / "Interpro_entry_mapping.pkl.gz"), result.stdout)
        self.assertEqual((interpro / "interpro.xml.gz").read_bytes(), self.xml_fixture.read_bytes())
        self.assertEqual((interpro / "interpro.xml.gz.sha256").read_text().strip(),
                         hashlib.sha256(self.xml_fixture.read_bytes()).hexdigest())

    def retained_interpro(self):
        xml = self.cache / "InterPro/interpro.xml.gz"
        xml.write_bytes(self.xml_fixture.read_bytes())
        Path(str(xml) + ".sha256").write_text(hashlib.sha256(xml.read_bytes()).hexdigest() + "\n")
        return xml

    def test_interpro_current_cache_reuses_exact_source_without_download_or_parser(self):
        xml = self.retained_interpro()
        mapping = self.cache / "InterPro/Interpro_entry_mapping.pkl.gz"
        mapping.write_bytes(b"existing mapping")
        newest = max(xml.stat().st_mtime, (ROOT / "scripts/protein_domain_mapping.py").stat().st_mtime) + 10
        os.utime(mapping, (newest, newest))
        result = self.run_shell('InterPro_parsing config')
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertNotIn("UNEXPECTED", result.stdout)
        self.assertNotIn("UPDATE:", result.stdout)
        self.assertEqual(mapping.read_bytes(), b"existing mapping")

    def test_interpro_missing_pickle_rebuild_uses_retained_source_without_download(self):
        xml = self.retained_interpro()
        before = xml.read_bytes()
        result = self.run_shell(r'''
            python() {
                printf 'PARSER:%s\n' "$*"
                while [[ $# -gt 0 ]]; do
                    if [[ "$1" == --mapping_output ]]; then printf 'rebuilt mapping\n' > "$2"; break; fi
                    shift
                done
            }
            InterPro_parsing config
        ''')
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertNotIn("UNEXPECTED_NETWORK", result.stdout)
        self.assertIn("PARSER:", result.stdout)
        self.assertEqual(xml.read_bytes(), before)

    def test_interpro_changed_retained_source_fails_without_replacing_previous_mapping(self):
        xml = self.retained_interpro()
        with gzip.open(xml, "wb") as stream:
            stream.write(b"<different/>\n")
        mapping = self.cache / "InterPro/Interpro_entry_mapping.pkl.gz"
        mapping.write_bytes(b"previous mapping")
        before = xml.read_bytes()
        result = self.run_shell('InterPro_parsing config')
        self.assertNotEqual(result.returncode, 0)
        self.assertIn("source identity differs", result.stderr)
        self.assertNotIn("UNEXPECTED", result.stdout)
        self.assertEqual(mapping.read_bytes(), b"previous mapping")
        self.assertEqual(xml.read_bytes(), before)

    def test_interpro_existing_invalid_or_empty_source_is_retained_without_refresh(self):
        xml = self.cache / "InterPro/interpro.xml.gz"
        for content in (b"not gzip", b""):
            with self.subTest(content=content):
                xml.write_bytes(content)
                result = self.run_shell('InterPro_parsing config')
                self.assertNotEqual(result.returncode, 0)
                self.assertNotIn("UNEXPECTED", result.stdout)
                self.assertNotIn("UPDATE:", result.stdout)
                self.assertEqual(xml.read_bytes(), content)

    def test_failed_clingen_erepo_download_stops_before_parser_and_dosage(self):
        result = self.run_shell(r'''
            wget() { printf 'DOWNLOAD_FAILED\n'; return 74; }
            ClinGen_deploy config
        ''')
        self.assertNotEqual(result.returncode, 0)
        self.assertEqual(result.stdout.count("DOWNLOAD_FAILED"), 1)
        self.assertNotIn("UNEXPECTED_PYTHON", result.stdout)
        self.assertNotIn("UPDATE:", result.stdout)

    def test_clingen_parser_failure_cannot_be_masked_by_dosage_success(self):
        self.clinvar_inputs()
        dosage = self.cache / "ClinGen/gene_dosage_sensitivity.hg19.tsv"
        before = dosage.read_bytes()
        result = self.run_shell(r'''
            python() { printf 'PARSER_FAILED\n'; return 75; }
            ClinGen_deploy config
        ''')
        self.assertNotEqual(result.returncode, 0)
        self.assertIn("PARSER_FAILED", result.stdout)
        self.assertNotIn("UNEXPECTED_NETWORK", result.stdout)
        self.assertNotIn("UPDATE:", result.stdout)
        self.assertEqual(dosage.read_bytes(), before)

    def test_clingen_map_config_failure_stops_before_dosage_processing(self):
        self.clinvar_inputs()
        result = self.run_shell(r'''
            python() { return 0; }
            update_yaml() { printf 'UPDATE:%s\n' "$2"; [[ "$2" != clingen_map ]]; }
            ClinGen_deploy config
        ''')
        self.assertNotEqual(result.returncode, 0)
        self.assertEqual(result.stdout.splitlines(), ["UPDATE:clingen_map"])
        self.assertNotIn("UNEXPECTED_NETWORK", result.stdout)

    def test_failed_or_empty_clingen_dosage_download_never_updates_dosage_config(self):
        self.clinvar_inputs(dosage=False)
        for download_status in (0, 76):
            with self.subTest(download_status=download_status):
                self.env["TEST_DOWNLOAD_STATUS"] = str(download_status)
                result = self.run_shell(r'''
                    python() { return 0; }
                    wget() { printf 'DOSAGE_DOWNLOAD\n'; return "$TEST_DOWNLOAD_STATUS"; }
                    ClinGen_deploy config
                ''')
                self.assertNotEqual(result.returncode, 0)
                self.assertIn("UPDATE:clingen_map ", result.stdout)
                self.assertNotIn("UPDATE:gene_dosage_sensitivity", result.stdout)

    def test_clingen_success_keeps_source_parser_filter_and_order(self):
        result = self.run_shell(r'''
            awk() { printf 'AWK:%s\n' "$*" >&2; command awk "$@"; }
            wget() {
                printf 'DOWNLOAD:%s\n' "$*"
                if [[ "${@: -1}" == *gene_dosage_sensitivity* ]]; then
                    printf '# comment\n#Gene Symbol\tGenomic Location\tHaploinsufficiency Score\tHaploinsufficiency Description\tHaploinsufficiency Disease ID\tTriplosensitivity Score\nGENE\tchr1:1-2\t3\tSufficient evidence\tMONDO:0000001\t0\n' > "${@: -1}"
                else
                    printf 'ClinVar Variation Id\tApplied Evidence Codes (Met)\n123\tPS3\n' > "${@: -1}"
                fi
            }
            python() {
                if [[ "$1" == *mechanism_source_input.py ]]; then "$TEST_PYTHON" "$@"; return; fi
                printf 'PARSER:%s\n' "$*"
                while [[ $# -gt 0 ]]; do
                    if [[ "$1" == -o ]]; then printf 'stub mapping\n' > "$2"; break; fi
                    shift
                done
            }
            ClinGen_deploy config
        ''')
        self.assertEqual(result.returncode, 0, result.stderr)
        clingen = self.cache / "ClinGen"
        self.assertIn("DOWNLOAD:-c https://erepo.clinicalgenome.org/evrepo/api/summary/classifications/download -O " +
                      str(clingen / "clingen_erepo.tsv.download."), result.stdout)
        self.assertIn("prepare_clingen_map.py -v " + str(self.folder / "clinvar.vcf.gz") +
                      " -t " + str(clingen / "clingen_erepo.tsv") +
                      " -o " + str(clingen / "clingen_map.hg19.pkl.gz"), result.stdout)
        self.assertIn("DOWNLOAD:-c http://ftp.clinicalgenome.org/ClinGen_gene_curation_list_GRCh37.tsv -O " +
                      str(clingen / "gene_dosage_sensitivity.hg19.tsv.download."), result.stdout)
        self.assertIn(r"AWK:-F \t $1 !~ /^#/ {print;} $1 ~ /^#Gene/{print;}", result.stderr)
        updates = [line.split(" ", 1)[0] for line in result.stdout.splitlines() if line.startswith("UPDATE:")]
        self.assertEqual(updates, ["UPDATE:clingen_map", "UPDATE:gene_dosage_sensitivity"])
        self.assertEqual((clingen / "gene_dosage_sensitivity.hg19.tsv").read_text(),
                         "#Gene Symbol\tGenomic Location\tHaploinsufficiency Score\tHaploinsufficiency Description\tHaploinsufficiency Disease ID\tTriplosensitivity Score\nGENE\tchr1:1-2\t3\tSufficient evidence\tMONDO:0000001\t0\n")


if __name__ == "__main__":
    unittest.main()

