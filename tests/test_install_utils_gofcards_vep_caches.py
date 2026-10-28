"""GoFCards prepares exact-release official caches only when its build needs them."""
import importlib.util
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile
import unittest

ROOT = Path(__file__).resolve().parents[1]
DATA_TMP = Path(os.environ.get("PRIVA_TEST_TMPDIR", tempfile.gettempdir())) / "priva-gofcards-cache-tests"
SPEC = importlib.util.spec_from_file_location("gofcards_prerequisites", ROOT / "scripts/installer_prerequisites.py")
PREREQS = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(PREREQS)


class GoFCardsMergedCaches(unittest.TestCase):
    def setUp(self):
        DATA_TMP.mkdir(parents=True, exist_ok=True)
        self.temp = tempfile.TemporaryDirectory(prefix=ROOT.name + "_", dir=DATA_TMP)
        self.addCleanup(self.temp.cleanup)
        self.folder = Path(self.temp.name)
        self.cache = self.folder / "cache"
        self.plugins = self.folder / "plugins"
        self.plugin_cache = self.folder / "plugin-cache"
        for path in (self.cache, self.plugins, self.plugin_cache):
            path.mkdir()
        self.references = {}
        for assembly, length in (("hg19", 249250621), ("hg38", 248956422)):
            path = self.folder / (assembly + ".fa")
            path.write_text(">chr1\nA\n")
            Path(str(path) + ".fai").write_text(f"chr1\t{length}\t6\t1\t2\n")
            self.references[assembly] = path
        self.chain = self.folder / "hg19ToHg38.over.chain"
        self.chain.write_text("chain 1 chr1 249250621 + 0 1 chr1 248956422 + 0 1 1\n1\n")
        self.hgnc = self.folder / "hgnc.tsv"
        self.hgnc.write_text("fixture\n")
        self.config = {"assembly": "hg38", "ref_genome": str(self.references["hg38"]),
                       "vep_installed_version": "116", "vep_cache_dir": str(self.cache),
                       "vep_plugins_dir": str(self.plugins), "vep_plugins_cachedir": str(self.plugin_cache),
                       "gofcards_hg19_fasta": str(self.references["hg19"]),
                       "gofcards_hg38_fasta": str(self.references["hg38"]),
                       "gofcards_hg19_to_hg38_chain": str(self.chain),
                       "gofcards_workdir": str(self.folder / "work"), "gofcards_exact_gof_cache": str(self.folder / "fresh_deployed.json.gz"), "hgnc_table": str(self.hgnc)}
        self.config_path = self.folder / "config.json"
        self.events = self.folder / "events.jsonl"
        self.api = self.folder / "fake-api.py"
        self.api.write_text('''import json,os,pathlib,sys
args=sys.argv[1:]
def option(name): return args[args.index(name)+1]
assembly=option('--VEP_ASSEMBLY')
assembly={'hg19':'GRCh37','hg38':'GRCh38'}.get(assembly,assembly)
release=option('--VEP_VERSION')
with open(os.environ['TEST_EVENTS'],'a') as f: f.write(json.dumps({'event':'CACHE','assembly':assembly,'release':release})+'\\n')
if os.environ.get('TEST_FAIL_ASSEMBLY')==assembly: raise SystemExit(73)
destination=pathlib.Path(option('--VEP_CACHEDIR'))/'homo_sapiens_merged'/(os.environ.get('TEST_WRONG_RELEASE',release)+'_'+assembly)
(destination/'1').mkdir(parents=True,exist_ok=True)
(destination/'info.txt').write_text('species\\thomo_sapiens\\nassembly\\t'+assembly+'\\nvar_type\\ttabix\\n')
for name in ('chr_synonyms.txt','1/all_vars.gz','1/all_vars.gz.csi'): (destination/name).write_text('fixture\\n')
''')
        self.env = dict(os.environ, TEST_PYTHON=sys.executable, TEST_CONFIG=str(self.config_path),
                        TEST_EVENTS=str(self.events), TEST_API=str(self.api), PYTHONDONTWRITEBYTECODE="1",
                        PRIVA_INSTALL_INSPECT_ROOTS=os.pathsep.join((str(ROOT), str(self.folder))))
        for key in ("PRIVA_SMOKE_SKIP_VEP", "PRIVA_FORCE_GOFCARDS_CACHE", "PRIVA_FORCE_ALL_CACHES"):
            self.env.pop(key, None)

    def seed(self, assembly, release="116", metadata_assembly=None):
        root = self.cache / "homo_sapiens_merged" / (release + "_" + assembly)
        (root / "1").mkdir(parents=True, exist_ok=True)
        (root / "info.txt").write_text("species\thomo_sapiens\nassembly\t" + (metadata_assembly or assembly) + "\nvar_type\ttabix\n")
        for name in ("chr_synonyms.txt", "1/all_vars.gz", "1/all_vars.gz.csi"):
            (root / name).write_text("fixture\n")
        return root

    def run_shell(self, body):
        self.config_path.write_text(json.dumps(self.config))
        prefix = r'''
source "$1"
log() { printf '%s\n' "$*" >&2; }
read_yaml() {
 "$TEST_PYTHON" -c 'import json,sys; v=json.load(open(sys.argv[1]));
for key in sys.argv[2].split("."): v=v.get(key) if isinstance(v,dict) else None
print("null" if v is None else v)' "$1" "$2"
}
update_yaml() { :; }
update_or_append_yaml() { :; }
priva_verify_vep_release() { [[ $1 == 116 ]]; }
vep_cache_api_install() { "$TEST_PYTHON" "$TEST_API" "$@"; }
vep() { :; }
bcftools() { :; }
CrossMap() { :; }
gofcards_reference_install() { printf '{"event":"REFERENCES"}\n' >> "$TEST_EVENTS"; }
validate_gofcards_exact_gof_cache() { return 0; }
wget() { printf 'UNEXPECTED_NETWORK\n'; return 97; }
curl() { printf 'UNEXPECTED_NETWORK\n'; return 97; }
python() {
 if [[ $1 == */build_gofcards_exact_gof_cache.py ]]; then
   printf '{"event":"BUILD"}\n' >> "$TEST_EVENTS"
   "$TEST_PYTHON" -c 'import json,sys; print("BUILD_ARGS="+json.dumps(sys.argv[1:]))' "$@"
 else "$TEST_PYTHON" "$@"; fi
}
'''
        return subprocess.run(["bash", "--noprofile", "--norc", "-c", prefix + body, "test",
                               str(ROOT / "scripts/install_utils.sh")], env=self.env,
                              text=True, capture_output=True, timeout=45)

    def recorded(self):
        return [json.loads(line) for line in self.events.read_text().splitlines()] if self.events.exists() else []

    def test_fresh_build_prepares_both_requested_release_caches_for_each_cohort_assembly(self):
        for assembly in ("hg19", "hg38"):
            with self.subTest(assembly=assembly):
                shutil.rmtree(self.cache)
                self.cache.mkdir()
                self.events.unlink(missing_ok=True)
                self.config.update(assembly=assembly, ref_genome=str(self.references[assembly]))
                result = self.run_shell('gofcards_exact_gof_cache_install "$TEST_CONFIG"')
                self.assertEqual(result.returncode, 0, result.stderr)
                self.assertEqual(self.recorded(), [{"event": "REFERENCES"},
                    {"event": "CACHE", "assembly": "GRCh37", "release": "116"},
                    {"event": "CACHE", "assembly": "GRCh38", "release": "116"}, {"event": "BUILD"}])
                arguments = json.loads(next(line.split("=", 1)[1] for line in result.stdout.splitlines() if line.startswith("BUILD_ARGS=")))
                self.assertEqual(arguments[arguments.index("--hg19-fasta") + 1], str(self.references["hg19"]))
                self.assertEqual(arguments[arguments.index("--hg38-fasta") + 1], str(self.references["hg38"]))
                self.assertEqual(arguments[arguments.index("--vep-cache-version") + 1], "116")
                self.assertEqual(json.loads(self.config_path.read_text()), self.config)

    def test_main_preflight_defers_only_absent_caches_without_creating_work_or_downloading(self):
        result = self.run_shell('gofcards_exact_gof_cache_install "$TEST_CONFIG" --main-preflight')
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertEqual(self.recorded(), [])
        self.assertFalse((self.folder / "work").exists())

    def test_complete_pair_is_reused_without_calling_cache_installer(self):
        for assembly in ("GRCh37", "GRCh38"):
            self.seed(assembly)
        result = self.run_shell('gofcards_vep_caches_install "$TEST_CONFIG"')
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertEqual(self.recorded(), [])
        self.assertEqual(json.loads(self.config_path.read_text()), self.config)

    def test_selected_complete_cache_only_acquires_missing_opposite_assembly(self):
        self.seed("GRCh38")
        result = self.run_shell('gofcards_vep_caches_install "$TEST_CONFIG"')
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertEqual(self.recorded(), [{"event": "CACHE", "assembly": "GRCh37", "release": "116"}])

    def test_current_normalized_artifact_reuses_without_any_cache_or_reference_preparation(self):
        target = self.folder / "work/gofcards_exact_gof.normalized.json.gz"
        target.parent.mkdir()
        target.write_text("validated by fixture\n")
        os.utime(target, (target.stat().st_atime + 300, target.stat().st_mtime + 300))
        self.references["hg19"].unlink()
        result = self.run_shell('gofcards_exact_gof_cache_install "$TEST_CONFIG"')
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertEqual(self.recorded(), [])
        self.assertIn("normalized cache is current", result.stderr)

    def test_partial_or_wrong_assembly_existing_cache_fails_before_any_download(self):
        for malformed in ("missing_index", "wrong_assembly"):
            with self.subTest(malformed=malformed):
                shutil.rmtree(self.cache)
                self.cache.mkdir()
                root = self.seed("GRCh38", metadata_assembly="GRCh37" if malformed == "wrong_assembly" else "GRCh38")
                if malformed == "missing_index":
                    (root / "1/all_vars.gz.csi").unlink()
                result = self.run_shell('gofcards_exact_gof_cache_install "$TEST_CONFIG" --main-preflight')
                self.assertNotEqual(result.returncode, 0)
                self.assertIn("GoFCards GRCh38 merged VEP 116 cache", result.stderr)
                self.assertEqual(self.recorded(), [])
                self.assertTrue(root.exists())

    def test_api_success_with_wrong_release_never_reaches_scientific_builder(self):
        self.env["TEST_WRONG_RELEASE"] = "115"
        result = self.run_shell('gofcards_exact_gof_cache_install "$TEST_CONFIG"')
        self.assertNotEqual(result.returncode, 0)
        self.assertNotIn({"event": "BUILD"}, self.recorded())
        self.assertIn("complete VEP 116 GRCh37", result.stderr)

    def test_explicit_existing_cache_parent_file_is_not_deferred(self):
        for location in ("cache", "merged_parent"):
            with self.subTest(location=location):
                if self.cache.is_dir():
                    shutil.rmtree(self.cache)
                else:
                    self.cache.unlink()
                if location == "cache":
                    destination = self.cache
                else:
                    self.cache.mkdir()
                    destination = self.cache / "homo_sapiens_merged"
                destination.write_text("retain malformed explicit input\n")
                result = self.run_shell('gofcards_exact_gof_cache_install "$TEST_CONFIG" --main-preflight')
                self.assertNotEqual(result.returncode, 0)
                self.assertIn("cache parent is not a directory", result.stderr)
                self.assertEqual(self.recorded(), [])
                self.assertEqual(destination.read_text(), "retain malformed explicit input\n")

    def test_other_installed_release_is_not_used_as_fallback(self):
        self.seed("GRCh37", "115")
        result = self.run_shell('gofcards_vep_caches_install "$TEST_CONFIG"')
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertEqual({entry["release"] for entry in self.recorded()}, {"116"})
        self.assertTrue((self.cache / "homo_sapiens_merged/115_GRCh37").is_dir())

    def test_smoke_mode_does_not_defer_or_acquire_missing_caches(self):
        self.env["PRIVA_SMOKE_SKIP_VEP"] = "1"
        result = self.run_shell('gofcards_exact_gof_cache_install "$TEST_CONFIG" --main-preflight')
        self.assertNotEqual(result.returncode, 0)
        self.assertEqual(self.recorded(), [])

    def test_standalone_validator_cannot_silently_defer_missing_caches(self):
        with self.assertRaisesRegex(ValueError, "merged VEP 116 cache"):
            PREREQS.validate_gofcards(self.config, main_preflight=True)

    def test_opposite_cache_failure_stops_main_install_before_success(self):
        self.env["TEST_FAIL_ASSEMBLY"] = "GRCh37"
        result = self.run_shell(r'''
conda_install_vep() { PRIVA_INSTALL_ENV_NAME=priva_acmg; }
priva_install_basic_paths() { PRIVA_INSTALL_TMPDIR="${TEST_CONFIG%/*}"; }
priva_manual_cache_inputs() { :; }
priva_prepare_install_config() { priva_install_prerequisites "$2"; }
priva_install_prerequisites() { gofcards_exact_gof_cache_install "$1" --main-preflight; }
bundled_snapshot_install() { :; }
nm_to_enst_mapping_install() { :; }
priva_configure_cache_paths() { :; }
for function in reference_genome_install VEP_plugins_install LoFtee_install gnomAD_install AlphaMissense_anno_gnomAD InterPro_parsing ClinVar_VCF_deploy ClinVar_AA_stat ClinVar_Gene_stat AlphaMissense_stat AlphaMissense_pick_intolerant_domains prepare_pm1_regions ClinVar_patho_AF_stat CADD_install CADD_runtime_install CADD_upstream_check ClinGen_deploy CDS_FASTA_install pext_install; do
 eval "$function() { :; }"
done
ClinVar_install() { [[ $# -eq 2 && $1 == --config && $2 == "$TEST_CONFIG" ]]; }
mechanism_resource_install() { gofcards_exact_gof_cache_install "$1"; }
main_install "$TEST_CONFIG"
''')
        self.assertNotEqual(result.returncode, 0)
        self.assertIn("Failed to install HPO/MONDO or mechanism resources", result.stderr)
        self.assertNotIn("installation is completed", result.stderr)
        self.assertNotIn({"event": "BUILD"}, self.recorded())
        calls = [entry for entry in self.recorded() if entry["event"] == "CACHE"]
        self.assertEqual(calls, [{"event": "CACHE", "assembly": "GRCh38", "release": "116"},
                                 {"event": "CACHE", "assembly": "GRCh37", "release": "116"}])


if __name__ == "__main__":
    unittest.main()
