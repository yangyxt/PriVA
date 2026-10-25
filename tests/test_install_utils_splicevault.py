"""SpliceVault deployment guards with local BGZF fixtures and no network."""
import gzip
import os
from pathlib import Path
import shutil
import subprocess
import tempfile
import unittest

ROOT = Path(__file__).resolve().parents[1]
INSTALLER = ROOT / "scripts/install_utils.sh"
TSV = (b"#chrom\tpos\tref\talt\tsplicevault_site_loc\n"
       b"chr1\t10\tA\tC\tchr1:11\nchr1\t20\tG\tT\tchr1:21\n")


@unittest.skipUnless(shutil.which("bgzip") and shutil.which("tabix"), "htslib tools required")
class SpliceVaultInstall(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory(prefix="priva_splicevault_",
                                               dir=os.environ.get("PRIVA_TEST_TMPDIR", os.environ.get("TMPDIR")))
        self.addCleanup(self.temp.cleanup)
        self.folder = Path(self.temp.name)
        self.cache = self.folder / "cache"
        self.reference = self.folder / "reference.fa"
        self.reference.write_text(">chr1\nA\n")
        self.chain = self.folder / "data/liftover/test.chain.gz"
        self.chain.parent.mkdir(parents=True)
        self.chain.write_text("stub chain\n")
        self.source = self.folder / "source.tsv.gz"
        with self.source.open("wb") as stream:
            subprocess.run(["bgzip", "-c"], input=TSV, stdout=stream, check=True)
        self.converted = self.folder / "converted.tsv"
        self.converted.write_bytes(
            b"#chrom\tpos\tref\talt\tsplicevault_site_loc\n"
            b"chr2\t20\tC\tT\tchr2:21\nchr1\t30\tA\tC\tchr1:31\n")
        self.env = dict(os.environ, TEST_FOLDER=str(self.folder), TEST_CACHE=str(self.cache),
                        TEST_REFERENCE=str(self.reference), TEST_CHAIN=str(self.chain),
                        TEST_SOURCE=str(self.source), TEST_CONVERTED=str(self.converted),
                        TEST_TMP=str(self.folder / "tmp/nested"),
                        TEST_PRESCORE=str(self.cache / "SpliceVault/SpliceVault_data_GRCh38.tsv.gz"),
                        TEST_URL="https://fixtures.invalid/SpliceVault_data_GRCh38.tsv.gz")

    def run_shell(self, body):
        settings = r'''
            DATA_DIR="$TEST_FOLDER/data"
            log() { printf '%s\n' "$*" >&2; }
            read_yaml() {
                case "$2" in
                    assembly) printf '%s\n' "${TEST_ASSEMBLY:-hg38}" ;;
                    splicevault_prescore) printf '%s\n' "$TEST_PRESCORE" ;;
                    splicevault_url) printf '%s\n' "$TEST_URL" ;;
                    hg38_hg19_chain) printf '%s\n' "$TEST_CHAIN" ;;
                    ref_genome) printf '%s\n' "$TEST_REFERENCE" ;;
                    tmp_dir) printf '%s\n' "$TEST_TMP" ;;
                esac
            }
            update_yaml() { printf 'UPDATE:%s\n' "$*"; }
            wget() { printf 'UNEXPECTED_NETWORK\n'; return 98; }
        '''
        return subprocess.run(["bash", "-c", 'source "$1"; ' + settings + body,
                               "test", str(INSTALLER)], text=True, capture_output=True, env=self.env)

    def local_download(self):
        return r'''
            wget() {
                local url destination
                printf 'DOWNLOAD:%s\n' "$*" >&2
                while [[ $# -gt 0 ]]; do
                    case "$1" in
                        -c) shift ;;
                        -O) destination=$2; shift 2 ;;
                        *) url=$1; shift ;;
                    esac
                done
                case "$url" in
                    "$TEST_URL") cp "$TEST_SOURCE" "$destination" ;;
                    "${TEST_URL}.tbi") cp "${TEST_SOURCE}.tbi" "$destination" ;;
                    *) return 99 ;;
                esac
            }
        '''

    def index_source(self):
        subprocess.run(["tabix", "-s", "1", "-b", "2", "-e", "2", "-c", "#",
                        str(self.source)], check=True)

    def test_existing_bgzf_without_index_is_indexed_before_reuse(self):
        self.env["TEST_PRESCORE"] = str(self.source)
        original = self.source.read_bytes()
        result = self.run_shell('SpliceVault_install config "$TEST_CACHE"')
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertNotIn("UNEXPECTED_NETWORK", result.stdout)
        self.assertNotIn("UPDATE", result.stdout)
        self.assertEqual(self.source.read_bytes(), original)
        self.assertTrue(Path(str(self.source) + ".tbi").is_file())
        records = subprocess.check_output(["tabix", str(self.source), "chr1:10-10"])
        self.assertEqual(records, b"chr1\t10\tA\tC\tchr1:11\n")

    def test_existing_current_index_is_reused_without_download(self):
        self.index_source()
        self.env["TEST_PRESCORE"] = str(self.source)
        index = Path(str(self.source) + ".tbi")
        before = (index.read_bytes(), index.stat().st_mtime_ns)
        result = self.run_shell('SpliceVault_install config "$TEST_CACHE"')
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertNotIn("UNEXPECTED_NETWORK", result.stdout)
        self.assertEqual((index.read_bytes(), index.stat().st_mtime_ns), before)

    def test_stale_or_corrupt_index_is_rebuilt(self):
        for corrupt in (False, True):
            with self.subTest(corrupt=corrupt):
                index = Path(str(self.source) + ".tbi")
                if index.exists():
                    index.unlink()
                self.index_source()
                if corrupt:
                    index.write_bytes(b"invalid index\n")
                else:
                    os.utime(index, (1, 1))
                self.env["TEST_PRESCORE"] = str(self.source)
                original = self.source.read_bytes()
                result = self.run_shell('SpliceVault_install config "$TEST_CACHE"')
                self.assertEqual(result.returncode, 0, result.stderr)
                self.assertNotIn("UNEXPECTED_NETWORK", result.stdout)
                self.assertEqual(self.source.read_bytes(), original)
                self.assertEqual(subprocess.check_output(["tabix", "-l", str(self.source)]), b"chr1\n")
                self.assertGreaterEqual(index.stat().st_mtime_ns, self.source.stat().st_mtime_ns)

    def test_invalid_existing_compression_fails_without_download(self):
        for ordinary_gzip in (False, True):
            with self.subTest(ordinary_gzip=ordinary_gzip):
                self.source.write_bytes(gzip.compress(TSV) if ordinary_gzip else b"invalid gzip\n")
                self.env["TEST_PRESCORE"] = str(self.source)
                result = self.run_shell('SpliceVault_install config "$TEST_CACHE"')
                self.assertNotEqual(result.returncode, 0)
                self.assertNotIn("UNEXPECTED_NETWORK", result.stdout)
                self.assertNotIn("UPDATE", result.stdout)

    def test_grch38_download_creates_resource_directory_and_validates_index(self):
        self.index_source()
        result = self.run_shell(self.local_download() + 'SpliceVault_install config "$TEST_CACHE"')
        self.assertEqual(result.returncode, 0, result.stderr)
        target = self.cache / "SpliceVault/SpliceVault_data_GRCh38.tsv.gz"
        self.assertEqual(target.read_bytes(), self.source.read_bytes())
        self.assertEqual(subprocess.check_output(["tabix", "-l", str(target)]), b"chr1\n")
        self.assertIn("UPDATE:config splicevault_prescore " + str(target), result.stdout)

    def test_corrupt_downloaded_index_never_updates_configuration(self):
        Path(str(self.source) + ".tbi").write_bytes(b"invalid index\n")
        result = self.run_shell(self.local_download() + 'SpliceVault_install config "$TEST_CACHE"')
        self.assertNotEqual(result.returncode, 0)
        self.assertNotIn("UPDATE", result.stdout)

    def test_grch19_conversion_keeps_original_command_and_sort_rules(self):
        self.env["TEST_ASSEMBLY"] = "hg19"
        self.env["TEST_PRESCORE"] = str(self.cache / "SpliceVault/SpliceVault_data_GRCh37.tsv.gz")
        result = self.run_shell(self.local_download() + r'''
            python() {
                if [[ "$1" != */SpliceVault_tsv_vcf_conversion.py ]]; then command python "$@"; return $?; fi
                printf 'CONVERT:%s\n' "$*" >&2
                local output
                while [[ $# -gt 0 ]]; do
                    if [[ "$1" == --output_tsv ]]; then output=$2; shift; fi
                    shift
                done
                cp "$TEST_CONVERTED" "$output"
            }
            sort() {
                printf 'SORT:%s\n' "$*" >&2
                local -a args=()
                while [[ $# -gt 0 ]]; do
                    [[ "$1" == --parallel=8 ]] || args+=("$1")
                    shift
                done
                command sort --parallel=1 "${args[@]}"
            }
            announce_remove_tmps() { rm -f "$@"; }
            SpliceVault_install config "$TEST_CACHE"
        ''')
        self.assertEqual(result.returncode, 0, result.stderr)
        target = self.cache / "SpliceVault/SpliceVault_data_GRCh37.tsv.gz"
        self.assertTrue(Path(self.env["TEST_TMP"]).is_dir())
        self.assertIn("SpliceVault_tsv_vcf_conversion.py", result.stderr)
        self.assertIn("--chain_file " + str(self.chain), result.stderr)
        self.assertIn("--reference_fasta " + str(self.reference), result.stderr)
        self.assertIn("--temp_dir " + self.env["TEST_TMP"], result.stderr)
        self.assertIn("SORT:-k1,1V -k2,2n --parallel=8 --buffer-size=8G -T " +
                      self.env["TEST_TMP"], result.stderr)
        self.assertEqual(gzip.decompress(target.read_bytes()),
                         b"#chrom\tpos\tref\talt\tsplicevault_site_loc\n"
                         b"chr1\t30\tA\tC\tchr1:31\nchr2\t20\tC\tT\tchr2:21\n")
        self.assertEqual(subprocess.check_output(["tabix", "-l", str(target)]), b"chr1\nchr2\n")
        self.assertIn("UPDATE:config splicevault_prescore " + str(target), result.stdout)

    def test_missing_or_empty_conversion_inputs_stop_before_download(self):
        self.env["TEST_ASSEMBLY"] = "hg19"
        self.env["TEST_PRESCORE"] = str(self.cache / "SpliceVault/SpliceVault_data_GRCh37.tsv.gz")
        for variable in ("TEST_REFERENCE", "TEST_CHAIN", "TEST_TMP"):
            with self.subTest(variable=variable):
                original = self.env[variable]
                if variable == "TEST_TMP":
                    self.env[variable] = ""
                else:
                    self.env[variable] = str(self.folder / "absent.input")
                result = self.run_shell('SpliceVault_install config "$TEST_CACHE"')
                self.assertNotEqual(result.returncode, 0)
                self.assertNotIn("UNEXPECTED_NETWORK", result.stdout)
                self.assertNotIn("UPDATE", result.stdout)
                self.env[variable] = original
        self.reference.write_bytes(b"")
        result = self.run_shell('SpliceVault_install config "$TEST_CACHE"')
        self.assertNotEqual(result.returncode, 0)
        self.assertNotIn("UNEXPECTED_NETWORK", result.stdout)

    def test_failed_conversion_stops_before_index_or_config_update(self):
        self.env["TEST_ASSEMBLY"] = "hg19"
        self.env["TEST_PRESCORE"] = str(self.cache / "SpliceVault/SpliceVault_data_GRCh37.tsv.gz")
        result = self.run_shell(self.local_download() + r'''
            python() { return 71; }
            tabix() { printf 'UNEXPECTED_TABIX\n'; }
            SpliceVault_install config "$TEST_CACHE"
        ''')
        self.assertNotEqual(result.returncode, 0)
        self.assertNotIn("UNEXPECTED_TABIX", result.stdout)
        self.assertNotIn("UPDATE", result.stdout)

    def test_missing_url_and_unsupported_assembly_stop_before_download(self):
        for key, value in (("TEST_URL", ""), ("TEST_ASSEMBLY", "GRCh39")):
            with self.subTest(key=key):
                original = self.env.get(key)
                self.env[key] = value
                result = self.run_shell('SpliceVault_install config "$TEST_CACHE"')
                self.assertNotEqual(result.returncode, 0)
                self.assertNotIn("UNEXPECTED_NETWORK", result.stdout)
                self.assertNotIn("UPDATE", result.stdout)
                if original is None:
                    self.env.pop(key)
                else:
                    self.env[key] = original

    def test_failed_download_never_updates_configuration(self):
        result = self.run_shell(r'''
            wget() { printf 'DOWNLOAD_FAILED\n'; return 72; }
            SpliceVault_install config "$TEST_CACHE"
        ''')
        self.assertNotEqual(result.returncode, 0)
        self.assertIn("DOWNLOAD_FAILED", result.stdout)
        self.assertNotIn("UPDATE", result.stdout)


if __name__ == "__main__":
    unittest.main()

