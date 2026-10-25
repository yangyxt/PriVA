"""No-acquisition fixtures for the existing config-filename whitespace contract.

Storage roots still reject spaces. Existing readable CONFIG filenames and the
selected reference remain single arguments across installer calls. Scientific
stage bodies are stubbed only in dispatcher tests; UTR reuse uses real YAML reads.
"""
import json
import os
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest
import yaml

ROOT = Path(__file__).resolve().parents[1]
INSTALLER = ROOT / "scripts/install_utils.sh"
STAGES = ["VEP_plugins_install", "LoFtee_install", "gnomAD_install", "AlphaMissense_anno_gnomAD",
          "InterPro_parsing", "pext_install", "ClinVar_install", "AlphaMissense_stat", "AlphaMissense_pick_intolerant_domains",
          "prepare_pm1_regions", "ClinVar_patho_AF_stat", "CADD_install", "CADD_runtime_install",
          "ClinGen_deploy", "CDS_FASTA_install", "mechanism_resource_install"]
PLUGINS = ["UTRAnnotator_install", "LOEUF_install", "Conservation_install", "AlphaMissense_install",
           "AlphaMissense_anno", "AlphaMissense_pick_missense_constrained_segments", "SpliceAI_install",
           "SpliceVault_install", "MaveDB_install", "PrimateAI_install"]
PRELUDE = r'''
source "$1"
wget() { echo UNEXPECTED_NETWORK >&2; return 98; }
curl() { echo UNEXPECTED_NETWORK >&2; return 98; }
record_call() {
    "$TEST_PYTHON" -c 'import os,json,sys;open(os.environ["CALLS"],"a").write(json.dumps(sys.argv[1:])+"\n")' "$@"
    [[ "$1" != "${FAIL_STAGE:-}" ]]
}
'''


class ConfigFilenameSpaces(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory(dir=os.environ.get("PRIVA_TEST_TMPDIR", os.environ.get("TMPDIR")))
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)
        self.ref = self.root / "selected reference.fa"
        self.ref.write_text(">chr1\nA\n")
        self.config = {"assembly": "hg38", "vep_installed_version": 116, "ref_genome": str(self.ref)}
        for key in ("cache_root", "tmp_dir", "vep_cache_dir", "vep_plugins_dir", "vep_plugins_cachedir"):
            path = self.root / key
            path.mkdir()
            self.config[key] = str(path)
        self.calls = self.root / "calls.jsonl"

    def execute(self, script, name="valid config.yaml", fail=""):
        config = self.root / name
        config.write_text(yaml.safe_dump(self.config, sort_keys=False, width=4096))
        self.calls.write_text("")
        env = {key: value for key, value in os.environ.items() if not key.startswith("PRIVA_SMOKE_")}
        env.update(TEST_PYTHON=sys.executable, CALLS=str(self.calls), CONFIG=str(config), FAIL_STAGE=fail,
                   CONDA_PREFIX=str(Path(sys.executable).parent.parent))
        result = subprocess.run(["bash", "-c", PRELUDE + script, "fixture", str(INSTALLER)],
                                env=env, text=True, capture_output=True, timeout=30)
        self.assertNotIn("UNEXPECTED_NETWORK", result.stdout + result.stderr)
        return result, config, [json.loads(line) for line in self.calls.read_text().splitlines()]

    def main_script(self):
        setup = r'''
conda_install_vep() { PRIVA_INSTALL_ENV_NAME=priva_acmg; }
priva_manual_cache_inputs() { :; }
priva_prepare_install_config() { :; }
bundled_snapshot_install() { :; }
reference_genome_install() { :; }
priva_verify_vep_release() { :; }
vep_cache_api_install() { :; }
'''
        return setup + "\n".join(name + '() { record_call ' + name + ' "$@"; }' for name in STAGES) + '\nmain_install "$CONFIG"\n'

    def test_main_dispatch_preserves_config_and_reference_for_spaces_and_plain_names(self):
        for filename in ("valid config.yaml", "ordinary.yaml"):
            with self.subTest(filename=filename):
                result, config, calls = self.execute(self.main_script(), filename)
                self.assertEqual(result.returncode, 0, result.stderr)
                self.assertEqual([row[0] for row in calls], STAGES)
                expected = {name: [str(config)] for name in STAGES}
                expected.update(VEP_plugins_install=[self.config['vep_plugins_dir'], self.config['vep_plugins_cachedir'], 'hg38', str(config), 'priva_acmg'],
                                LoFtee_install=[self.config['vep_plugins_cachedir'], str(config)],
                                gnomAD_install=[str(config), str(self.ref)], ClinVar_install=['--config', str(config)],
                                CADD_install=[str(config), 'hg38'], CADD_runtime_install=[str(config), 'both'])
                for name, *args in calls:
                    self.assertEqual(args, expected[name], name)

    def test_main_stage_failure_stops_following_stages_with_exact_config(self):
        result, config, calls = self.execute(self.main_script(), fail="gnomAD_install")
        self.assertNotEqual(result.returncode, 0)
        self.assertEqual([row[0] for row in calls], STAGES[:3])
        self.assertEqual(calls[-1][1:], [str(config), str(self.ref)])

    def test_pext_ready_before_clinvar_consumer_fixture(self):
        # Actual main dispatch; acquisition and ClinVar/VEP bodies are stubs.
        # This marker models track readiness, not scientific annotation.
        overrides = r'''
pext_install() {
    record_call pext_install "$@" || return $?
    printf 'fixture installed track\n' > "$CALLS.pext-ready"
}
ClinVar_install() {
    [[ -s "$CALLS.pext-ready" ]] || return 81
    record_call ClinVar_install "$@"
}
'''
        script = self.main_script().replace('\nmain_install "$CONFIG"\n', overrides + '\nmain_install "$CONFIG"\n')
        result, _, calls = self.execute(script)
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertEqual([row[0] for row in calls], STAGES)
        self.assertEqual(sum(row[0] == "pext_install" for row in calls), 1)

    def test_pext_failure_stops_before_clinvar_consumer_fixture(self):
        result, _, calls = self.execute(self.main_script(), fail="pext_install")
        self.assertNotEqual(result.returncode, 0)
        self.assertEqual([row[0] for row in calls], STAGES[:STAGES.index("pext_install") + 1])
        self.assertNotIn("ClinVar_install", [row[0] for row in calls])

    def test_plugin_dispatch_preserves_config_and_argument_order(self):
        script = '\n'.join(name + '() { record_call ' + name + ' "$@"; }' for name in PLUGINS)
        script += '\nVEP_plugins_install "' + self.config['vep_plugins_dir'] + '" "' + self.config['vep_plugins_cachedir'] + '" hg38 "$CONFIG" priva_acmg\n'
        result, config, calls = self.execute(script)
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertEqual([row[0] for row in calls], PLUGINS)
        cache, plugins = self.config['vep_plugins_cachedir'], self.config['vep_plugins_dir']
        expected = [[str(config), cache], [str(config), 'hg38'], [str(config), cache, 'hg38'],
                    [str(config), cache, 'hg38', plugins], [str(config)], [str(config)],
                    [str(config), cache, plugins], [str(config), cache, str(ROOT/'scripts')],
                    [str(config), cache], [str(config), cache, plugins]]
        self.assertEqual([row[1:] for row in calls], expected)
        result, config, calls = self.execute(script, fail='UTRAnnotator_install')
        self.assertNotEqual(result.returncode, 0)
        self.assertEqual(len(calls), 1)
        self.assertEqual(calls[0][1:], [str(config), cache])

    def test_real_utr_reads_selected_input_from_spaced_config(self):
        selected = self.root / 'selected utr.tsv'
        selected.write_text('opaque selected nonempty fixture\n')
        self.config['utr_annotator_file'] = str(selected)
        before = selected.read_bytes(), selected.stat().st_mtime_ns
        result, config, calls = self.execute('UTRAnnotator_install "$CONFIG" "' + self.config['vep_plugins_cachedir'] + '"\n')
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertEqual(yaml.safe_load(config.read_text()), self.config)
        self.assertEqual((selected.read_bytes(), selected.stat().st_mtime_ns), before)
        self.assertEqual(calls, [])

    def test_storage_root_whitespace_is_still_rejected_before_dispatch(self):
        path = self.root / 'unsupported cache root'
        path.mkdir()
        self.config['cache_root'] = str(path)
        result, _, calls = self.execute(self.main_script())
        self.assertNotEqual(result.returncode, 0)
        self.assertEqual(calls, [])
        self.assertIn('without spaces or shell characters', result.stderr)


if __name__ == '__main__':
    unittest.main()
