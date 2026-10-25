"""Bundled LOEUF needs existence and nonzero size only, for the selected assembly."""
import os
from pathlib import Path
import subprocess
import tempfile
import unittest


ROOT = Path(__file__).resolve().parents[1]


class BundledLoeuf(unittest.TestCase):
    def setUp(self):
        temp = tempfile.TemporaryDirectory(prefix="priva_bundled_loeuf_", dir=os.environ.get("PRIVA_TEST_TMPDIR", os.environ.get("TMPDIR")))
        self.addCleanup(temp.cleanup)
        self.folder = Path(temp.name)
        self.data = self.folder / "source data's" / "data"
        (self.data / "loeuf").mkdir(parents=True)
        self.config = self.folder / "selected config's.yaml"
        self.config.write_text("loeuf_prescore: ''\n")

    def run_loeuf(self, assembly):
        body = r'''
source "$1" >/dev/null || exit
DATA_DIR="$TEST_DATA"
update_yaml() {
    [[ $# == 3 && $1 == "$TEST_CONFIG" ]] || exit 97
    printf 'REGISTER %s\n' "$3"
}
gzip() { printf 'UNEXPECTED_CONTENT_CHECK\n'; exit 98; }
python() { printf 'UNEXPECTED_DECODE\n'; exit 99; }
LOEUF_install "$TEST_CONFIG" "$TEST_ASSEMBLY"
'''
        return subprocess.run(["bash", "--noprofile", "--norc", "-c", body, "test", str(ROOT / "scripts/install_utils.sh")],
                              env=dict(os.environ, TEST_DATA=str(self.data), TEST_CONFIG=str(self.config), TEST_ASSEMBLY=assembly),
                              text=True, capture_output=True, timeout=15)

    def test_nonzero_opaque_bundle_registers_selected_assembly_and_alias(self):
        for name in ("loeuf_dataset.tsv.gz", "loeuf_dataset_hg38.tsv.gz"):
            (self.data / "loeuf" / name).write_bytes(b"opaque bundled data\n")
        for assembly, name in (("hg19", "loeuf_dataset.tsv.gz"), ("GRCh37", "loeuf_dataset.tsv.gz"),
                               ("hg38", "loeuf_dataset_hg38.tsv.gz"), ("GRCh38", "loeuf_dataset_hg38.tsv.gz")):
            with self.subTest(assembly=assembly):
                result = self.run_loeuf(assembly)
                self.assertEqual(result.returncode, 0, result.stderr)
                self.assertEqual(result.stdout.strip(), "REGISTER " + str(self.data / "loeuf" / name))

    def test_missing_selected_bundle_does_not_register_or_substitute_other_assembly(self):
        (self.data / "loeuf/loeuf_dataset_hg38.tsv.gz").write_bytes(b"other assembly bundle\n")
        result = self.run_loeuf("hg19")
        self.assertNotEqual(result.returncode, 0)
        self.assertIn("missing or empty", result.stderr)
        self.assertNotIn("REGISTER", result.stdout)

    def test_explicit_nonempty_selection_is_preserved_instead_of_bundle(self):
        (self.data / "loeuf/loeuf_dataset_hg38.tsv.gz").write_bytes(b"bundled snapshot")
        selected = self.folder / "selected.tsv.gz"
        selected.write_bytes(b"explicit opaque selected bytes")
        self.config.write_text("loeuf_prescore: " + str(selected) + "\n")
        result = self.run_loeuf("hg38")
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertEqual(result.stdout.strip(), "REGISTER " + str(selected))
        self.assertEqual(selected.read_bytes(), b"explicit opaque selected bytes")

    def test_missing_explicit_selection_never_falls_back_to_bundle(self):
        (self.data / "loeuf/loeuf_dataset_hg38.tsv.gz").write_bytes(b"bundled snapshot")
        self.config.write_text("loeuf_prescore: " + str(self.folder / "absent.tsv.gz") + "\n")
        result = self.run_loeuf("hg38")
        self.assertNotEqual(result.returncode, 0)
        self.assertNotIn("REGISTER", result.stdout)
        self.assertIn("no replacement selected", result.stderr)

    def test_empty_selected_bundle_does_not_register(self):
        (self.data / "loeuf/loeuf_dataset.tsv.gz").touch()
        result = self.run_loeuf("GRCh37")
        self.assertNotEqual(result.returncode, 0)
        self.assertIn("missing or empty", result.stderr)
        self.assertNotIn("REGISTER", result.stdout)


if __name__ == "__main__":
    unittest.main()
