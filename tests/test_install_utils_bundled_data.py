"""Exact canonical bundles reuse opaque nonzero bytes; external builds retain their routes."""
import importlib.util
import json
import os
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest
from unittest import mock


ROOT = Path(__file__).resolve().parents[1]
SPEC = importlib.util.spec_from_file_location("bundled_portable_paths", ROOT / "scripts/installer_config_paths.py")
PATHS = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(PATHS)

PRELUDE = r'''
source "$1" >/dev/null || exit
BASE_DIR="$TEST_SOURCE"
DATA_DIR="$BASE_DIR/data"
read_yaml() {
    "$TEST_PYTHON" -c 'import json,sys; v=json.load(open(sys.argv[1])).get(sys.argv[2]); print("null" if v is None else v)' "$1" "$2"
}
update_yaml() { printf 'REGISTER %s %s\n' "$2" "$3"; }
update_or_append_yaml() { update_yaml "$@"; }
python() { printf 'UNEXPECTED_BUILD %s\n' "$*"; return 91; }
wget() { printf 'UNEXPECTED_DOWNLOAD %s\n' "$*"; return 92; }
priva_verify_vep_release() { printf 'UNEXPECTED_VEP_CHECK\n'; return 93; }
'''


class BundledDataPolicy(unittest.TestCase):
    def setUp(self):
        temp = tempfile.TemporaryDirectory(prefix="priva_bundled_policy_", dir=os.environ.get("PRIVA_TEST_TMPDIR", os.environ.get("TMPDIR")))
        self.addCleanup(temp.cleanup)
        self.folder = Path(temp.name)
        self.source = self.folder / "checkout"
        self.cache = self.folder / "managed_cache"
        self.config_path = self.folder / "config.json"
        self.config = {"assembly": "hg19", "vep_installed_version": "116",
                       "cache_root": str(self.cache), "tmp_dir": str(self.folder),
                       "gene_mechanism_cache_dir": str(self.cache / "patho_mechanism"),
                       "gene_mechanism_raw_dir": str(self.cache / "patho_mechanism/raw"),
                       "gofcards_workdir": str(self.cache / "gofcards/build_work"),
                       "clinvar_vcv_xml": str(self.cache / "patho_mechanism/raw/clinvar_vcv/ClinVarVCVRelease_00-latest_weekly.xml.gz")}
        self.bundles = {}
        for key, relative in PATHS.BUNDLED_CACHES.items():
            path = self.source / relative
            path.parent.mkdir(parents=True, exist_ok=True)
            path.write_bytes(b"opaque nonempty original bundle, not decoded by installer\n")
            self.bundles[key] = path
            self.config[key] = str(path)
        for assembly in ("hg19", "hg38"):
            for key, name in (("clingen_map", f"clingen_map.{assembly}.pkl.gz"),
                              ("gene_dosage_sensitivity", f"gene_dosage_sensitivity.{assembly}.tsv")):
                path = self.source / "data/clingen" / name
                path.parent.mkdir(parents=True, exist_ok=True)
                path.write_bytes(b"opaque original ClinGen bundle\n")
                self.bundles[key + "." + assembly] = path
                if assembly == "hg19":
                    self.config[key] = str(path)
        erepo = self.source / "data/clingen/clingen_erepo.tsv"
        erepo.write_bytes(b"original bundled eRepo source\n")
        self.bundles["clingen_erepo"] = erepo
        self.env = {key: value for key, value in os.environ.items() if not key.startswith("PRIVA_FORCE_")}
        self.env.update(TEST_SOURCE=str(self.source), TEST_CONFIG=str(self.config_path),
                        TEST_PYTHON=sys.executable, TEST_FOLDER=str(self.folder), PYTHONDONTWRITEBYTECODE="1")

    def run_shell(self, body, **env):
        self.config_path.write_text(json.dumps(self.config))
        return subprocess.run(["bash", "--noprofile", "--norc", "-c", PRELUDE + body,
                               "test", str(ROOT / "scripts/install_utils.sh")],
                              text=True, capture_output=True, env=dict(self.env, **env), timeout=30)

    def snapshot(self):
        return {key: path.read_bytes() for key, path in self.bundles.items()}

    def test_all_six_canonical_routes_reuse_opaque_bundles_without_build_inputs(self):
        before = self.snapshot()
        for command in ('gene_pathogenic_mechanism_cache_install "$TEST_CONFIG"',
                        'gofcards_exact_gof_cache_install "$TEST_CONFIG" --main-preflight',
                        'gofcards_exact_gof_cache_install "$TEST_CONFIG"',
                        'gofcards_clinvar_injection_install "$TEST_CONFIG"',
                        'gene_nonlof_mechanism_cache_install "$TEST_CONFIG"',
                        'hpo_condition_mechanism_cache_install "$TEST_CONFIG"'):
            with self.subTest(command=command):
                result = self.run_shell(command)
                self.assertEqual(result.returncode, 0, result.stderr)
                self.assertNotIn("UNEXPECTED", result.stdout)
        self.assertEqual(self.snapshot(), before)
        self.assertFalse(self.cache.exists())

    def test_sources_only_skips_only_when_both_final_bundles_are_selected(self):
        result = self.run_shell('gene_nonlof_mechanism_cache_install "$TEST_CONFIG" sources-only')
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertNotIn("UNEXPECTED", result.stdout)
        self.assertFalse(self.cache.exists())

    def test_external_gof_output_still_fetches_vcv_with_bundled_nonlof(self):
        self.config["gofcards_exact_gof_cache"] = str(self.cache / "gofcards/custom.json.gz")
        builder = self.folder / "source_only_builder.py"
        builder.write_text("# source-only command fixture\n")
        schema = self.folder / "schema.json"
        schema.write_text("{}\n")
        hgnc = self.folder / "hgnc.tsv"
        hgnc.write_text("opaque source fixture\n")
        self.config.update(gene_nonlof_mechanism_builder_script=str(builder),
                           gene_nonlof_mechanism_schema=str(schema), hgnc_table=str(hgnc))
        result = self.run_shell(r'''
            python() {
                printf 'FETCH_ARGS %s\n' "$*"
                "$TEST_PYTHON" - "$TEST_CONFIG" <<'PY'
import json,sys
from pathlib import Path
path=Path(json.load(open(sys.argv[1]))['clinvar_vcv_xml'])
path.parent.mkdir(parents=True,exist_ok=True)
path.write_bytes(b'operational source-fetch fixture')
PY
            }
            gene_nonlof_mechanism_cache_install "$TEST_CONFIG" sources-only
        ''')
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertIn("--fetch-sources-only", result.stdout)
        self.assertIn("FETCH_ARGS " + str(builder), result.stdout)
        self.assertTrue(Path(self.config["clinvar_vcv_xml"]).is_file())

    def test_mechanism_integration_keeps_stage_order_without_rebuilding_bundles(self):
        result = self.run_shell('''
            mondo_hpo_scope_install() { printf 'SCOPE_STAGE\n'; }
            mechanism_resource_install "$TEST_CONFIG"
        ''')
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertEqual(result.stdout.strip(), "SCOPE_STAGE")
        text = (ROOT / "scripts/install_utils.sh").read_text()
        block = text.split("function mechanism_resource_install() {", 1)[1].split("\n}", 1)[0]
        stages = ["mondo_hpo_scope_install", "gene_pathogenic_mechanism_cache_install",
                  "gofcards_exact_gof_cache_install", "gene_nonlof_mechanism_cache_install",
                  "gofcards_clinvar_injection_install", "gene_nonlof_mechanism_cache_install",
                  "hpo_condition_mechanism_cache_install"]
        previous = -1
        for stage in stages:
            previous = block.index("    " + stage + " ", previous + 1)

    def test_force_cannot_overwrite_explicit_canonical_bundle(self):
        before = self.snapshot()
        for function in ("gene_pathogenic_mechanism_cache_install", "gofcards_exact_gof_cache_install",
                         "gofcards_clinvar_injection_install", "gene_nonlof_mechanism_cache_install",
                         "hpo_condition_mechanism_cache_install", "ClinGen_deploy"):
            with self.subTest(function=function):
                result = self.run_shell(function + ' "$TEST_CONFIG"', PRIVA_FORCE_ALL_CACHES="1")
                self.assertNotEqual(result.returncode, 0)
                self.assertIn("explicit data-side output", result.stderr)
                self.assertNotIn("UNEXPECTED", result.stdout)
        self.assertEqual(self.snapshot(), before)
        self.assertFalse(self.cache.exists())

    def test_unset_outputs_select_bundles_and_force_selects_managed_outputs(self):
        config = {key: self.config[key] for key in ("cache_root", "tmp_dir", "assembly", "vep_installed_version")}
        force_env = {key: "0" for key in os.environ if key.startswith("PRIVA_FORCE_")}
        with mock.patch.dict(os.environ, force_env):
            updates = PATHS.portable_updates(config, self.source)
        for key in PATHS.BUNDLED_CACHES:
            self.assertEqual(updates[key], str(self.bundles[key]))
        self.assertEqual(updates["clingen_map"], str(self.bundles["clingen_map.hg19"]))
        with mock.patch.dict(os.environ, PRIVA_FORCE_ALL_CACHES="1"):
            updates = PATHS.portable_updates(config, self.source)
        for key in (*PATHS.BUNDLED_CACHES, "clingen_map", "gene_dosage_sensitivity"):
            self.assertTrue(Path(updates[key]).is_relative_to(self.cache), key)

    def test_missing_empty_or_external_link_bundle_defaults_to_managed_output(self):
        config = {key: self.config[key] for key in ("cache_root", "tmp_dir", "assembly", "vep_installed_version")}
        for state in ("missing", "empty", "external_link"):
            with self.subTest(state=state):
                selected = self.bundles["gene_nonlof_mechanism_json"]
                if selected.exists() or selected.is_symlink():
                    selected.unlink()
                if state == "empty":
                    selected.touch()
                if state == "external_link":
                    external = self.folder / "external_nonlof.json.gz"
                    external.write_bytes(b"external bytes retain their policy")
                    selected.symlink_to(external)
                updates = PATHS.portable_updates(config, self.source)
                self.assertEqual(updates["gene_nonlof_mechanism_json"], str(self.cache / "patho_mechanism/gene_nonlof_mechanism_curated_assertions.json.gz"))

    def test_empty_canonical_file_is_not_reused_or_decoded(self):
        selected = self.bundles["ddg2p_mechanism_evidence"]
        selected.write_bytes(b"")
        result = self.run_shell('priva_bundled_cache_reuse --key ddg2p_mechanism_evidence --path "' + str(selected) + '"')
        self.assertEqual(result.returncode, 1)
        self.assertNotIn("UNEXPECTED", result.stdout)
        self.assertEqual(selected.read_bytes(), b"")

    def test_external_symlink_cannot_gain_bundle_trust_and_keeps_validation(self):
        selected = self.bundles["ddg2p_mechanism_evidence"]
        selected.unlink()
        external = self.folder / "external_evidence.tsv"
        external.write_bytes(b"external original bytes\n")
        selected.symlink_to(external)
        self.config["gene_mechanism_builder_script"] = str(ROOT / "scripts/build_gene_pathogenic_mechanism_cache.py")
        result = self.run_shell('''
            python() { printf 'EXTERNAL_BUILD\n'; }
            validate_gene_pathogenic_mechanism_evidence() { printf 'EXTERNAL_VALIDATE\n'; return 87; }
            gene_pathogenic_mechanism_cache_install "$TEST_CONFIG"
        ''')
        self.assertNotEqual(result.returncode, 0)
        self.assertIn("EXTERNAL_BUILD", result.stdout)
        self.assertIn("EXTERNAL_VALIDATE", result.stdout)
        self.assertEqual(external.read_bytes(), b"external original bytes\n")

    def test_clingen_both_bundles_reuse_for_both_assembly_aliases_without_parsing(self):
        before = self.snapshot()
        for assembly, canonical in (("hg19", "hg19"), ("GRCh37", "hg19"), ("hg38", "hg38"), ("GRCh38", "hg38")):
            with self.subTest(assembly=assembly):
                self.config.update(assembly=assembly, clingen_map=str(self.bundles["clingen_map." + canonical]),
                                   gene_dosage_sensitivity=str(self.bundles["gene_dosage_sensitivity." + canonical]))
                result = self.run_shell('ClinGen_deploy "$TEST_CONFIG"')
                self.assertEqual(result.returncode, 0, result.stderr)
                self.assertNotIn("UNEXPECTED", result.stdout)
        self.assertEqual(self.snapshot(), before)

    def test_clingen_bundled_map_is_untouched_when_only_external_dosage_is_built(self):
        before = self.snapshot()
        dosage = self.cache / "ClinGen/external_dosage.tsv"
        self.config["gene_dosage_sensitivity"] = str(dosage)
        result = self.run_shell(r'''
            wget() { printf '#Gene Symbol\tGenomic Location\tHaploinsufficiency Score\tHaploinsufficiency Description\tHaploinsufficiency Disease ID\tTriplosensitivity Score\nGENE1\tchr1:1-2\t3\tSufficient evidence\tMONDO:0000001\t0\n' > "${@: -1}"; }
            python() {
                [[ $1 == *mechanism_source_input.py ]] || { printf 'UNEXPECTED_BUILD %s\n' "$*"; return 91; }
                "$TEST_PYTHON" "$@"
            }
            ClinGen_deploy "$TEST_CONFIG"
        ''')
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertNotIn("UNEXPECTED_BUILD", result.stdout)
        self.assertTrue(dosage.is_file())
        self.assertEqual(self.snapshot(), before)

    def test_clingen_external_map_uses_bundled_erepo_without_touching_bundled_dosage(self):
        before = self.snapshot()
        selected = self.cache / "ClinGen/external_map.pkl.gz"
        self.config["clingen_map"] = str(selected)
        os.utime(self.bundles["clingen_erepo"], (1, 1))
        result = self.run_shell(r'''
            python() {
                printf 'MAP_ARGS %s\n' "$*"
                while [[ $# -gt 0 ]]; do
                    if [[ $1 == -o ]]; then printf 'new map fixture\n' > "$2"; return; fi
                    shift
                done
                return 99
            }
            ClinGen_deploy "$TEST_CONFIG"
        ''')
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertIn("-t " + str(self.bundles["clingen_erepo"]), result.stdout)
        self.assertNotIn("UNEXPECTED_DOWNLOAD", result.stdout)
        self.assertTrue(selected.is_file())
        self.assertEqual(self.snapshot(), before)

    def test_existing_external_erepo_keeps_precedence_over_bundled_source(self):
        selected = self.cache / "ClinGen/external_map.pkl.gz"
        selected.parent.mkdir(parents=True)
        erepo = selected.parent / "clingen_erepo.tsv"
        erepo.write_bytes(b"existing external eRepo identity\n")
        self.config["clingen_map"] = str(selected)
        result = self.run_shell(r'''
            python() {
                printf 'MAP_ARGS %s\n' "$*"
                while [[ $# -gt 0 ]]; do
                    if [[ $1 == -o ]]; then printf 'new map fixture\n' > "$2"; return; fi
                    shift
                done
                return 99
            }
            ClinGen_deploy "$TEST_CONFIG"
        ''')
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertIn("-t " + str(erepo), result.stdout)
        self.assertNotIn("-t " + str(self.bundles["clingen_erepo"]), result.stdout)
        self.assertNotIn("UNEXPECTED_DOWNLOAD", result.stdout)
        self.assertEqual(erepo.read_bytes(), b"existing external eRepo identity\n")


    def test_dated_canonical_erepo_reuses_when_map_is_missing_or_externally_linked(self):
        selected = self.bundles["clingen_map.hg19"]
        erepo = self.bundles["clingen_erepo"]
        original = erepo.read_bytes()
        os.utime(erepo, (1, 1))
        prior = erepo.stat()
        for state in ("missing", "external_link"):
            with self.subTest(state=state):
                if selected.exists() or selected.is_symlink():
                    selected.unlink()
                external = self.folder / "supplied_external_map.pkl.gz"
                if state == "external_link":
                    external.write_bytes(b"external map retains its own validation route\n")
                    selected.symlink_to(external)
                result = self.run_shell(r'''
                    python() { printf 'MAP_ARGS %s\n' "$*"; }
                    ClinGen_deploy "$TEST_CONFIG"
                ''')
                self.assertEqual(result.returncode, 0, result.stderr)
                self.assertIn("-t " + str(erepo), result.stdout)
                self.assertNotIn("UNEXPECTED_DOWNLOAD", result.stdout)
                self.assertEqual(erepo.read_bytes(), original)
                after = erepo.stat()
                self.assertEqual((after.st_ino, after.st_size, after.st_mtime_ns),
                                 (prior.st_ino, prior.st_size, prior.st_mtime_ns))
                if state == "external_link":
                    self.assertEqual(external.read_bytes(), b"external map retains its own validation route\n")


if __name__ == "__main__":
    unittest.main()
