"""Selected plugin inputs stay selected; acquisition is limited to managed defaults."""
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile
import unittest

ROOT = Path(__file__).resolve().parents[1]
INSTALLER = ROOT / "scripts/install_utils.sh"
PRELUDE = r'''
source "$1"
log() { printf '%s\n' "$*" >&2; }
read_yaml() {
    "$TEST_PYTHON" -c 'import json,sys; v=json.load(open(sys.argv[1])).get(sys.argv[2], ""); print(str(v).lower() if isinstance(v,bool) else v)' "$1" "$2"
}
update_yaml() {
    printf 'UPDATE:%s\n' "$2"
    "$TEST_PYTHON" -c 'import json,sys; p,k,v=sys.argv[1:]; d=json.load(open(p)); d[k]=v; open(p,"w").write(json.dumps(d))' "$@"
}
python() { "$TEST_PYTHON" "$@"; }
wget() { printf 'UNEXPECTED_NETWORK\n'; return 98; }
git() { printf 'UNEXPECTED_GIT\n'; return 97; }
'''


@unittest.skipUnless(shutil.which("bgzip") and shutil.which("tabix"), "htslib tools required")
class ExplicitPluginInputs(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory(
            prefix="priva_selected_plugin_", dir=os.environ.get("PRIVA_TEST_TMPDIR", os.environ.get("TMPDIR")))
        self.addCleanup(self.temporary.cleanup)
        self.folder = Path(self.temporary.name)
        self.cache = self.folder / "managed"
        self.cache.mkdir()
        self.config_file = self.folder / "config.json"
        self.config = {"assembly": "hg38", "vep_plugins_cachedir": str(self.cache)}
        self.env = dict(os.environ, TEST_PYTHON=sys.executable, TEST_CONFIG=str(self.config_file),
                        TEST_CACHE=str(self.cache))

    def run_shell(self, command, setup=""):
        self.config_file.write_text(json.dumps(self.config))
        return subprocess.run(["bash", "-c", PRELUDE + setup + "\n" + command, "test", str(INSTALLER)],
                              env=self.env, cwd=self.folder, text=True, capture_output=True, timeout=45)

    def bgzf(self, content, name="selected.tsv.gz"):
        path = self.folder / name
        with path.open("wb") as output:
            subprocess.run(["bgzip", "-c"], input=content, stdout=output, check=True)
        return path

    def test_mavedb_selected_table_repairs_only_missing_index_without_network(self):
        source = self.bgzf(b"#chrom\tstart\tend\nchr1\t1\t1\nchr1\t2\t2\n")
        self.config["mavedb_file"] = str(source)
        before = source.read_bytes(), source.stat().st_ctime_ns
        result = self.run_shell('MaveDB_install "$TEST_CONFIG" "$TEST_CACHE"')
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertNotIn("UNEXPECTED", result.stdout)
        self.assertNotIn("UPDATE", result.stdout)
        self.assertEqual((source.read_bytes(), source.stat().st_ctime_ns), before)
        self.assertEqual(json.loads(self.config_file.read_text()), self.config)
        self.assertEqual(subprocess.check_output(["tabix", str(source), "chr1:1-1"]), b"chr1\t1\t1\n")

    def test_missing_custom_table_inputs_fail_before_network_or_configuration_change(self):
        for key, command in (
            ("utr_annotator_file", 'UTRAnnotator_install "$TEST_CONFIG" "$TEST_CACHE"'),
            ("splicevault_prescore", 'SpliceVault_install "$TEST_CONFIG" "$TEST_CACHE"'),
            ("mavedb_file", 'MaveDB_install "$TEST_CONFIG" "$TEST_CACHE"'),
            ("primateai_prescore", 'PrimateAI_install "$TEST_CONFIG" "$TEST_CACHE" plugins'),
        ):
            with self.subTest(key=key):
                self.config[key] = str(self.cache / "arbitrary-missing.tsv.gz")
                result = self.run_shell(command)
                self.assertNotEqual(result.returncode, 0)
                self.assertNotIn("UNEXPECTED", result.stdout)
                self.assertNotIn("UPDATE", result.stdout)
                self.assertEqual(json.loads(self.config_file.read_text()), self.config)
                self.config.pop(key)

    def test_primate_selected_prepared_relative_path_reuses_and_repairs_index(self):
        source = self.bgzf(b"#chrom\tpos\tref\talt\tscore\nchr1\t1\tA\tC\t0.1\n")
        self.config["primateai_prescore"] = source.name
        before = source.read_bytes(), source.stat().st_ctime_ns
        result = self.run_shell('PrimateAI_install "$TEST_CONFIG" "$TEST_CACHE" plugins')
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertNotIn("UPDATE", result.stdout)
        self.assertEqual((source.read_bytes(), source.stat().st_ctime_ns), before)
        self.assertEqual(json.loads(self.config_file.read_text()), self.config)

    def test_zero_byte_selected_inputs_are_not_treated_as_ready(self):
        source = self.folder / "empty.tsv.gz"
        source.touch()
        for key, command in (
            ("primateai_prescore", 'PrimateAI_install "$TEST_CONFIG" "$TEST_CACHE" plugins'),
            ("mavedb_file", 'MaveDB_install "$TEST_CONFIG" "$TEST_CACHE"'),
        ):
            with self.subTest(key=key):
                self.config[key] = str(source)
                result = self.run_shell(command)
                self.assertNotEqual(result.returncode, 0)
                self.assertNotIn("UNEXPECTED", result.stdout)
                self.assertEqual(source.read_bytes(), b"")

    def test_blank_selections_register_existing_managed_defaults_without_network(self):
        for key, relative, command, content in (
            ("utr_annotator_file", "UTRannotator/uORF_5UTR_GRCh38_PUBLIC.txt",
             'UTRAnnotator_install "$TEST_CONFIG" "$TEST_CACHE"', b"opaque existing UTR resource"),
            ("mavedb_file", "MaveDB/MaveDB_variants.hg38.tsv.gz",
             'MaveDB_install "$TEST_CONFIG" "$TEST_CACHE"', None),
            ("splicevault_prescore", "SpliceVault/SpliceVault_data_GRCh38.tsv.gz",
             'SpliceVault_install "$TEST_CONFIG" "$TEST_CACHE"', None),
        ):
            with self.subTest(key=key):
                target = self.cache / relative
                target.parent.mkdir(parents=True)
                if content is not None:
                    target.write_bytes(content)
                else:
                    target.write_bytes(self.bgzf(b"#chrom\tstart\tend\nchr1\t1\t1\n").read_bytes())
                self.config[key] = ""
                result = self.run_shell(command)
                self.assertEqual(result.returncode, 0, result.stderr)
                self.assertNotIn("UNEXPECTED", result.stdout)
                self.assertEqual(json.loads(self.config_file.read_text())[key], str(target))
                self.config.pop(key)

    def test_mavedb_download_publishes_only_after_indexed_staging_succeeds(self):
        source = self.bgzf(b"#chrom\tstart\tend\nchr1\t1\t1\nchr1\t2\t2\n")
        subprocess.run(["tabix", "-s", "1", "-b", "2", "-e", "3", str(source)], check=True)
        self.env["TEST_DOWNLOAD_SOURCE"] = str(source)
        self.config["mavedb_url"] = "https://fixtures.invalid/mavedb.tsv.gz"
        setup = r'''
wget() {
    local url destination
    while [[ $# -gt 0 ]]; do
        case "$1" in -c) shift ;; -O) destination=$2; shift 2 ;; *) url=$1; shift ;; esac
    done
    [[ "$destination" == */.download.*/* ]] || return 99
    if [[ "$url" == *.tbi ]]; then cp "${TEST_DOWNLOAD_SOURCE}.tbi" "$destination"
    else cp "$TEST_DOWNLOAD_SOURCE" "$destination"; fi
}
'''
        result = self.run_shell('MaveDB_install "$TEST_CONFIG" "$TEST_CACHE"', setup)
        self.assertEqual(result.returncode, 0, result.stderr)
        target = self.cache / "MaveDB/MaveDB_variants.hg38.tsv.gz"
        self.assertEqual(target.read_bytes(), source.read_bytes())
        self.assertEqual(json.loads(self.config_file.read_text())["mavedb_file"], str(target))
        self.assertEqual(subprocess.check_output(["tabix", str(target), "chr1:2-2"]), b"chr1\t2\t2\n")

    def test_mavedb_failed_index_download_does_not_publish_data_or_register(self):
        setup = r'''
wget() {
    local destination=${@: -1}
    [[ "$destination" != *.tbi ]] || return 72
    printf 'partial fresh download' > "$destination"
}
'''
        self.config["mavedb_url"] = "https://fixtures.invalid/mavedb.tsv.gz"
        result = self.run_shell('MaveDB_install "$TEST_CONFIG" "$TEST_CACHE"', setup)
        self.assertNotEqual(result.returncode, 0)
        self.assertNotIn("UPDATE", result.stdout)
        self.assertFalse((self.cache / "MaveDB/MaveDB_variants.hg38.tsv.gz").exists())

    def test_interrupted_utr_clone_keeps_final_path_free_and_retry_succeeds(self):
        partial = r'''
git() { mkdir -p UTRannotator; printf 'interrupted clone' > UTRannotator/partial; return 71; }
'''
        result = self.run_shell('UTRAnnotator_install "$TEST_CONFIG" "$TEST_CACHE"', partial)
        self.assertNotEqual(result.returncode, 0)
        self.assertFalse((self.cache / "UTRannotator").exists())
        stages = list(self.cache.glob(".UTRannotator.clone.*"))
        self.assertEqual(len(stages), 1)
        self.assertEqual((stages[0] / "UTRannotator/partial").read_text(), "interrupted clone")
        complete = r'''
git() {
    mkdir -p UTRannotator
    printf 'publisher GRCh38 table' > UTRannotator/uORF_5UTR_GRCh38_PUBLIC.txt
    printf 'publisher GRCh37 table' > UTRannotator/uORF_5UTR_GRCh37_PUBLIC.txt
}
'''
        result = self.run_shell('UTRAnnotator_install "$TEST_CONFIG" "$TEST_CACHE"', complete)
        self.assertEqual(result.returncode, 0, result.stderr)
        target = self.cache / "UTRannotator/uORF_5UTR_GRCh38_PUBLIC.txt"
        self.assertEqual(target.read_text(), "publisher GRCh38 table")
        self.assertEqual(json.loads(self.config_file.read_text())["utr_annotator_file"], str(target))
        self.assertEqual(list(self.cache.glob(".UTRannotator.clone.*")), stages)

    def prepare_splicevault_liftover(self, rows, *, reference_name="reference.fa",
                                     chain_name="ten_bases.chain"):
        """Use actual converter/CrossMap and HTS tools; only acquisition is forbidden."""
        import pysam
        reference = self.folder / reference_name
        reference.write_text(">chr1\n" + "A" * 100 + "\n")
        pysam.faidx(str(reference))
        chain = self.folder / chain_name
        chain.write_text("chain 1 chr1 100 + 0 10 chr1 100 + 0 10 1\n10\n\n")
        directory = self.cache / "SpliceVault"
        directory.mkdir()
        self.splice_header = "#chrom\tpos\tref\talt\tsplicevault_site_loc\n"
        self.bgzf((self.splice_header + rows).encode(),
                  name=str(directory / "SpliceVault_data_GRCh38.tsv.gz"))
        self.config.update(assembly="hg19", splicevault_prescore="",
                           splicevault_url="https://invalid.example/no-network",
                           ref_genome=str(reference), hg38_hg19_chain=str(chain),
                           tmp_dir=str(self.folder / "liftover-tmp"))
        return directory / "SpliceVault_data_GRCh37.tsv.gz"

    @unittest.skipUnless(shutil.which("CrossMap"), "real CrossMap required")
    def test_splicevault_cold_unmapped_only_is_retained_but_not_registered(self):
        import gzip
        target = self.prepare_splicevault_liftover("chr1\t50\tA\tC\tchr1:50\n")
        result = self.run_shell('SpliceVault_install "$TEST_CONFIG" "$TEST_CACHE"')
        self.assertNotEqual(result.returncode, 0, result.stderr)
        self.assertNotIn("UPDATE", result.stdout)
        self.assertNotIn("UNEXPECTED", result.stdout)
        self.assertEqual(json.loads(self.config_file.read_text()), self.config)
        self.assertEqual(gzip.open(target, "rt").read(), self.splice_header)
        before = target.read_bytes(), target.stat().st_mtime_ns
        warm = self.run_shell('SpliceVault_install "$TEST_CONFIG" "$TEST_CACHE"')
        self.assertNotEqual(warm.returncode, 0, warm.stderr)
        self.assertNotIn("UPDATE", warm.stdout)
        self.assertEqual((target.read_bytes(), target.stat().st_mtime_ns), before)
        self.assertEqual(json.loads(self.config_file.read_text()), self.config)

    @unittest.skipUnless(shutil.which("CrossMap"), "real CrossMap required")
    def test_splicevault_real_mapped_cold_and_warm_preserve_sorted_rows(self):
        import gzip
        target = self.prepare_splicevault_liftover(
            "chr1\t8\tA\tC\tchr1:5\nchr1\t50\tA\tC\tchr1:50\nchr1\t3\tA\tG\tchr1:4\n")
        result = self.run_shell('SpliceVault_install "$TEST_CONFIG" "$TEST_CACHE"')
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertNotIn("UNEXPECTED", result.stdout)
        expected = self.splice_header + "chr1\t3\tA\tG\tchr1:4\nchr1\t8\tA\tC\tchr1:5\n"
        self.assertEqual(gzip.open(target, "rt").read(), expected)
        query = subprocess.check_output(["tabix", str(target), "chr1:1-10"])
        self.assertEqual(query.decode(), expected[len(self.splice_header):])
        self.config = json.loads(self.config_file.read_text())
        self.assertEqual(self.config["splicevault_prescore"], str(target))
        before = target.read_bytes(), target.stat().st_mtime_ns
        warm = self.run_shell('SpliceVault_install "$TEST_CONFIG" "$TEST_CACHE"')
        self.assertEqual(warm.returncode, 0, warm.stderr)
        self.assertNotIn("UPDATE", warm.stdout)
        self.assertNotIn("UNEXPECTED", warm.stdout)
        self.assertEqual((target.read_bytes(), target.stat().st_mtime_ns), before)
        self.assertEqual(json.loads(self.config_file.read_text()), self.config)
        self.assertEqual(subprocess.check_output(["tabix", str(target), "chr1:1-10"]), query)

    @unittest.skipUnless(shutil.which("CrossMap"), "real CrossMap required")
    def test_splicevault_real_selected_reference_and_chain_with_spaces_and_quotes(self):
        """Actual converter/CrossMap/HTS routing; no full-reference qualification."""
        import gzip
        target = self.prepare_splicevault_liftover(
            "chr1\t8\tA\tC\tchr1:5\nchr1\t50\tA\tC\tchr1:50\nchr1\t3\tA\tG\tchr1:4\n",
            reference_name="selected reference's.fa", chain_name="selected chain's.chain")
        sources = [Path(self.config[key]) for key in ("ref_genome", "hg38_hg19_chain")]
        original = {path: (path.read_bytes(), path.stat().st_ino, path.stat().st_mtime_ns)
                    for path in sources}
        result = self.run_shell('SpliceVault_install "$TEST_CONFIG" "$TEST_CACHE"')
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertNotIn("UNEXPECTED", result.stdout)
        expected = self.splice_header + "chr1\t3\tA\tG\tchr1:4\nchr1\t8\tA\tC\tchr1:5\n"
        self.assertEqual(gzip.open(target, "rt").read(), expected)
        query = subprocess.check_output(["tabix", str(target), "chr1:1-10"])
        self.assertEqual(query.decode(), expected[len(self.splice_header):])
        self.config = json.loads(self.config_file.read_text())
        self.assertEqual(self.config["splicevault_prescore"], str(target))
        before = target.read_bytes(), target.stat().st_ino, target.stat().st_mtime_ns
        warm = self.run_shell('SpliceVault_install "$TEST_CONFIG" "$TEST_CACHE"')
        self.assertEqual(warm.returncode, 0, warm.stderr)
        self.assertNotIn("UPDATE", warm.stdout)
        self.assertNotIn("UNEXPECTED", warm.stdout)
        self.assertEqual((target.read_bytes(), target.stat().st_ino, target.stat().st_mtime_ns), before)
        self.assertEqual(subprocess.check_output(["tabix", str(target), "chr1:1-10"]), query)
        self.assertEqual({path: (path.read_bytes(), path.stat().st_ino, path.stat().st_mtime_ns)
                          for path in sources}, original)

    def test_pext_validation_and_parser_preserve_caller_variables(self):
        command = r'''check_scope() {
            local tissues="Cells_EBV_transformedlymphocytes,Spleen,WholeBlood" tissue="caller-tissue"
            local expected=$tissues
            validate_pext_tissues "$tissues" hg38 || return 1
            [[ "$tissues" == "$expected" && "$tissue" == caller-tissue ]] || return 2
            get_pext_tissues_array "$tissues" > "$TEST_CACHE/parsed.txt" || return 3
            [[ "$tissues" == "$expected" && "$tissue" == caller-tissue ]] || return 4
            validate_pext_tissues "" hg38 || return 5
            get_pext_tissues_array "" > "$TEST_CACHE/empty.txt" || return 6
            [[ "$tissues" == "$expected" && "$tissue" == caller-tissue ]] || return 7
        }
        check_scope'''
        result = self.run_shell(command)
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertEqual((self.cache / "parsed.txt").read_text().splitlines(),
                         ["Cells_EBV_transformedlymphocytes", "Spleen", "WholeBlood"])
        self.assertEqual((self.cache / "empty.txt").read_bytes(), b"")

    def test_pext_empty_selection_acquires_only_default_mean(self):
        self.config.update(pext_enabled=True, pext_tissues="", pext_mean_bw="", pext_tissue_bw_dir="")
        setup = r'''wget() { [[ "$1" == -c && "$3" == -O ]] || return 9;
            printf fixture-mean > "$4"; }'''
        result = self.run_shell('pext_install "$TEST_CONFIG"', setup=setup)
        self.assertEqual(result.returncode, 0, result.stderr)
        directory = self.cache / "pext/hg38"
        self.assertEqual((directory / "exp_prop_mean.bw").read_bytes(), b"fixture-mean")
        self.assertFalse((directory / "tissues").exists())
        self.assertEqual(json.loads(self.config_file.read_text())["pext_tissues"], "")

    def test_pext_cold_route_acquires_every_selected_tissue_in_order(self):
        self.config.update(pext_enabled=True, pext_tissues="Cells_EBV_transformedlymphocytes, Spleen,WholeBlood",
                           pext_mean_bw="", pext_tissue_bw_dir="")
        setup = r'''wget() {
            printf 'DOWNLOAD:%s\n' "$2"
            local output
            while [[ $# -gt 0 ]]; do
                if [[ "$1" == -O ]]; then output=$2; shift 2; else shift; fi
            done
            printf fixture-resource > "$output"
        }'''
        result = self.run_shell('pext_install "$TEST_CONFIG"', setup=setup)
        self.assertEqual(result.returncode, 0, result.stderr)
        names = ["exp_prop_mean.bw", "Cells_EBV_transformedlymphocytes.bw", "Spleen.bw", "WholeBlood.bw"]
        requests = [line.removeprefix("DOWNLOAD:") for line in result.stdout.splitlines() if line.startswith("DOWNLOAD:")]
        self.assertEqual(requests, ["https://hgdownload.soe.ucsc.edu/gbdb/hg38/gnomAD/pext/" + name for name in names])
        directory = self.cache / "pext/hg38"
        for index, name in enumerate(names):
            self.assertEqual((directory / ("tissues" if index else "") / name).read_bytes(), b"fixture-resource")
        revised = json.loads(self.config_file.read_text())
        self.assertEqual(revised["pext_tissues"], self.config["pext_tissues"])
        self.assertEqual(revised["pext_mean_bw"], str(directory / names[0]))

    def test_pext_reuses_configured_mean_and_tissue_paths_without_acquisition(self):
        mean = self.folder / "custom mean.bw"
        mean.write_bytes(b"selected mean fixture")
        tissue_root = self.folder / "custom tissues"
        (tissue_root / "tissues").mkdir(parents=True)
        (tissue_root / "tissues/Lung.bw").write_bytes(b"selected Lung fixture")
        self.config.update(pext_enabled=True, pext_tissues="Lung", pext_mean_bw=str(mean),
                           pext_tissue_bw_dir=str(tissue_root))
        result = self.run_shell('pext_install "$TEST_CONFIG"')
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertNotIn("UNEXPECTED", result.stdout)
        self.assertNotIn("UPDATE", result.stdout)
        self.assertEqual(json.loads(self.config_file.read_text()), self.config)
        self.assertFalse((self.cache / "pext").exists())

    def test_pext_missing_custom_tissue_stops_before_default_mean_download(self):
        root = self.folder / "custom tissues"
        root.mkdir()
        self.config.update(pext_enabled=True, pext_tissues="Lung", pext_tissue_bw_dir=str(root))
        result = self.run_shell('pext_install "$TEST_CONFIG"')
        self.assertNotEqual(result.returncode, 0)
        self.assertNotIn("UNEXPECTED", result.stdout)
        self.assertNotIn("UPDATE", result.stdout)
        self.assertFalse((self.cache / "pext").exists())

    def test_invalid_pext_tissue_request_stops_before_any_acquisition(self):
        self.config.update(pext_enabled=True, pext_tissues="invalid-tissue")
        result = self.run_shell('pext_install "$TEST_CONFIG"')
        self.assertNotEqual(result.returncode, 0)
        self.assertNotIn("UNEXPECTED", result.stdout)
        self.assertNotIn("UPDATE", result.stdout)
        self.assertFalse((self.cache / "pext").exists())

    def test_disabled_pext_does_not_consume_missing_selected_inputs(self):
        self.config.update(pext_enabled=False, pext_mean_bw="missing", pext_tissue_bw_dir="missing")
        result = self.run_shell('pext_install "$TEST_CONFIG"')
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertNotIn("UNEXPECTED", result.stdout)
        self.assertNotIn("UPDATE", result.stdout)


if __name__ == "__main__":
    unittest.main()
