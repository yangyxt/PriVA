"""Conservation shell barriers with real BigWigs and a stubbed network writer."""
import json
import os
from pathlib import Path
import re
import subprocess
import sys
import tempfile
import unittest

from test_conservation_resource_check import guard, make_bigwig, make_reference


HERE = Path(__file__).resolve().parent
if (HERE / "Conservation_install.sh").is_file():
    FUNCTION = (HERE / "Conservation_install.sh").read_text()
    SCRIPT_DIR = HERE
else:
    SCRIPT_DIR = HERE.parent / "scripts"
    source = (SCRIPT_DIR / "install_utils.sh").read_text()
    FUNCTION = re.search(r"^function Conservation_install\(\).*?^}", source, re.M | re.S).group()


PRELUDE = r'''
set -euo pipefail
log() { printf '%s\n' "$*" >&2; }
event() { printf '%s\n' "$*" >> "$TEST_EVENTS"; }
read_yaml() {
    "$TEST_PYTHON" -c 'import json,sys; value=json.load(open(sys.argv[1])).get(sys.argv[2]); print("null" if value is None else value)' "$1" "$2"
}
python() { "$TEST_PYTHON" "$@"; }
wget() {
    event "wget $2"
    local output=${@: -1}
    cp -- "$TEST_DOWNLOAD" "$output"
    [[ "$TEST_FAILURE" != download ]] || return 71
    if [[ "$TEST_FAILURE" == publication ]]; then printf 'retained racing target' > "$TEST_TARGET"; fi
}
update_yaml() {
    event "yaml $2 $3"
    [[ "$TEST_FAILURE" != yaml ]]
}
'''


class ConservationShell(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory(
            prefix="priva_conservation_shell_", dir=os.environ.get("PRIVA_TEST_TMPDIR", os.environ.get("TMPDIR")))
        self.addCleanup(self.temporary.cleanup)
        self.folder = Path(self.temporary.name)
        self.cache = self.folder / "cache"
        self.target = self.cache / "Conservation/gerp_conservation_scores.homo_sapiens.GRCh38.bw"
        self.reference = make_reference(self.folder, "hg38")
        self.track = make_bigwig(self.folder, "valid_download", {"chr1": 248956422})
        self.config = {"assembly": "hg38", "conservation_file": str(self.target),
                       "ref_genome": str(self.reference),
                       "conservation_GRCh38_url": "https://example.invalid/explicit-custom.bw"}
        self.config_file = self.folder / "config.json"
        self.events = self.folder / "events.txt"
        self.events.touch()
        self.environment = dict(os.environ, SCRIPT_DIR=str(SCRIPT_DIR), TEST_CONFIG=str(self.config_file),
                                TEST_EVENTS=str(self.events), TEST_PYTHON=sys.executable,
                                TEST_CACHE=str(self.cache), TEST_DOWNLOAD=str(self.track),
                                TEST_TARGET=str(self.target), TEST_FAILURE="none", TEST_ASSEMBLY="hg38")

    def run_installer(self, failure="none"):
        self.config_file.write_text(json.dumps(self.config))
        environment = dict(self.environment, TEST_FAILURE=failure)
        result = subprocess.run(["bash", "-c", PRELUDE + FUNCTION + '\nConservation_install "$TEST_CONFIG" "$TEST_CACHE" "$TEST_ASSEMBLY"'],
                                env=environment, text=True, capture_output=True, timeout=30)
        return result, self.events.read_text().splitlines()

    def test_configured_custom_file_reuses_own_identity_without_default_release_claim(self):
        self.config["conservation_file"] = str(self.track)
        self.config["conservation_GRCh38_url"] = guard.GRCH38_DEFAULT_URL
        before = self.track.read_bytes()
        result, events = self.run_installer()
        self.assertEqual(result.returncode, 0, result.stderr)
        record = json.loads(result.stdout)
        self.assertEqual(events, [])
        self.assertIsNone(record["release"])
        self.assertIsNone(record["source_url"])
        self.assertEqual(self.track.read_bytes(), before)

    def test_invalid_configured_file_blocks_network_and_yaml_changes(self):
        self.track.write_bytes(b"retained corrupt configured file")
        self.config["conservation_file"] = str(self.track)
        result, events = self.run_installer()
        self.assertNotEqual(result.returncode, 0)
        self.assertEqual(events, [])
        self.assertEqual(self.track.read_bytes(), b"retained corrupt configured file")

    def test_missing_explicit_custom_file_is_not_replaced_by_a_default(self):
        self.config["conservation_file"] = str(self.folder / "missing-custom.bw")
        result, events = self.run_installer()
        self.assertNotEqual(result.returncode, 0)
        self.assertIn("Missing explicit Conservation input", result.stderr)
        self.assertEqual(events, [])
        self.assertFalse(self.target.exists())

    def test_fresh_explicit_download_validates_and_publishes_before_yaml(self):
        result, events = self.run_installer()
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertEqual(events, ["wget " + self.config["conservation_GRCh38_url"],
                                  "yaml conservation_file " + str(self.target)])
        self.assertEqual(self.target.read_bytes(), self.track.read_bytes())
        self.assertFalse(Path(str(self.target) + ".partial").exists())

    def test_default_download_requires_publisher_pin_before_publication_or_yaml(self):
        self.config["conservation_GRCh38_url"] = None
        result, events = self.run_installer()
        self.assertNotEqual(result.returncode, 0)
        self.assertEqual(events, ["wget " + guard.GRCH38_DEFAULT_URL])
        self.assertIn("publisher byte count mismatch", result.stderr)
        self.assertTrue(Path(str(self.target) + ".partial").is_file())
        self.assertFalse(self.target.exists())

    def test_failed_network_preserves_partial_and_blocks_yaml(self):
        result, events = self.run_installer("download")
        self.assertNotEqual(result.returncode, 0)
        self.assertEqual(len(events), 1)
        self.assertTrue(Path(str(self.target) + ".partial").is_file())
        self.assertFalse(self.target.exists())

    def test_corrupt_download_retains_partial_without_final_or_yaml(self):
        self.track.write_text("successful HTTP200 text body\n")
        result, events = self.run_installer()
        self.assertNotEqual(result.returncode, 0)
        self.assertEqual(len(events), 1)
        self.assertEqual(Path(str(self.target) + ".partial").read_bytes(), self.track.read_bytes())
        self.assertFalse(self.target.exists())

    def test_wrong_assembly_download_blocks_publication(self):
        wrong = make_bigwig(self.folder, "wrong", {"chr1": 249250621})
        self.environment["TEST_DOWNLOAD"] = str(wrong)
        result, events = self.run_installer()
        self.assertNotEqual(result.returncode, 0)
        self.assertEqual(len(events), 1)
        self.assertIn("chr1 length does not identify GRCh38", result.stderr)
        self.assertFalse(self.target.exists())
        self.assertTrue(Path(str(self.target) + ".partial").is_file())

    def test_publication_race_preserves_target_and_partial_without_yaml(self):
        result, events = self.run_installer("publication")
        self.assertNotEqual(result.returncode, 0)
        self.assertEqual(len(events), 1)
        self.assertEqual(self.target.read_bytes(), b"retained racing target")
        self.assertEqual(Path(str(self.target) + ".partial").read_bytes(), self.track.read_bytes())

    def test_yaml_update_failure_does_not_report_completed_deployment(self):
        result, events = self.run_installer("yaml")
        self.assertNotEqual(result.returncode, 0)
        self.assertEqual(len(events), 2)
        self.assertTrue(self.target.is_file())
        self.assertNotIn("is deployed to", result.stderr)

    def test_existing_canonical_target_validates_before_registering_without_download(self):
        self.target.parent.mkdir(parents=True)
        self.target.write_bytes(self.track.read_bytes())
        self.config["conservation_file"] = ""
        result, events = self.run_installer()
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertEqual(events, ["yaml conservation_file " + str(self.target)])
        self.assertIsNone(json.loads(result.stdout)["release"])

    def test_hg19_download_source_is_unchanged_and_native_validation_runs(self):
        self.target = self.cache / "Conservation/gerp_conservation_scores.homo_sapiens.GRCh37.bw"
        self.reference = make_reference(self.folder, "hg19")
        track = make_bigwig(self.folder, "hg19", {"chr1": 249250621})
        self.config.update(conservation_file=str(self.target), ref_genome=str(self.reference))
        self.environment.update(TEST_DOWNLOAD=str(track), TEST_TARGET=str(self.target), TEST_ASSEMBLY="GRCh37")
        result, events = self.run_installer()
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertEqual(events[0], "wget https://hgdownload.soe.ucsc.edu/gbdb/hg19/bbi/All_hg19_RS.bw")
        self.assertEqual(json.loads(result.stdout)["assembly"], "GRCh37")


if __name__ == "__main__":
    unittest.main()
