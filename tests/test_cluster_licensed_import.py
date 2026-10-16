"""Preflight barriers must run before reading any licensed cache file."""
import json
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest

ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / 'scripts/import_cluster_licensed.sbatch'
BODY = SCRIPT.read_text().split("<<'PY'\n", 1)[1].split('\nPY\n', 1)[0]


class LicensedImportPreflight(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory(prefix='licensed_import_', dir=ROOT / 'tests')
        self.addCleanup(self.temporary.cleanup)
        self.folder = Path(self.temporary.name)
        self.template = self.folder / 'template.yaml'
        self.template.write_text('assembly: hg19\n')
        self.config = self.folder / 'separate.yaml'
        self.confirmation = self.folder / 'confirmed.json'

    def run_import(self, confirmation):
        self.confirmation.write_text(json.dumps(confirmation))
        return subprocess.run(
            [sys.executable, '-c', BODY, str(ROOT), '054b4944055054148005016192647a8ed0384eef',
             str(self.template), str(self.config), str(self.confirmation), str(self.folder)],
            text=True, capture_output=True)

    def test_incomplete_acquisition_stops_before_template_or_cache_access(self):
        self.template.unlink()
        result = self.run_import({'status': 'downloading', 'confirmed_by': 'owner'})
        self.assertNotEqual(result.returncode, 0)
        self.assertIn('Explicit acquisition-owner completion confirmation', result.stderr)
        self.assertNotIn('FileNotFoundError', result.stderr)
        self.assertFalse(self.config.exists())

    def test_missing_owner_confirmation_is_rejected(self):
        result = self.run_import({'status': 'downloads_complete'})
        self.assertNotEqual(result.returncode, 0)
        self.assertIn('Explicit acquisition-owner completion confirmation', result.stderr)
        self.assertFalse(self.config.exists())

    def test_wrong_assembly_stops_before_cache_access(self):
        result = self.run_import({'status': 'downloads_complete', 'confirmed_by': 'owner',
                                  'assembly': 'hg38', 'files': {}})
        self.assertNotEqual(result.returncode, 0)
        self.assertIn('Confirmed files belong to another assembly', result.stderr)
        self.assertFalse(self.config.exists())

    def test_missing_resource_roles_stop_before_cache_access(self):
        result = self.run_import({'status': 'downloads_complete', 'confirmed_by': 'owner',
                                  'assembly': 'GRCh37', 'files': {}})
        self.assertNotEqual(result.returncode, 0)
        self.assertIn('All three matching-assembly resources are required', result.stderr)
        self.assertFalse(self.config.exists())

    def test_existing_config_is_preserved(self):
        self.config.write_text('existing configuration\n')
        result = self.run_import({'status': 'downloads_complete', 'confirmed_by': 'owner',
                                  'assembly': 'hg19', 'files': {}})
        self.assertNotEqual(result.returncode, 0)
        self.assertIn('Refusing to replace an existing import config', result.stderr)
        self.assertEqual(self.config.read_text(), 'existing configuration\n')


if __name__ == '__main__':
    unittest.main()
