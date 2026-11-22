"""Real installer/state helper, recording fixture priva and injected failures.

These tests establish orchestration/recovery, not domain model acceptance.
"""
import hashlib
import json
import os
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest
from types import SimpleNamespace
from unittest.mock import patch
import importlib.util

import fake_priva_binary


ROOT = Path(__file__).resolve().parents[1]
_spec = importlib.util.spec_from_file_location('das_install_state_test', ROOT / 'scripts/das_install_state.py')
STATE = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(STATE)


def builder(name):
    """The identity the fixture priva records for NAME."""
    return hashlib.sha256(name.encode()).hexdigest()


class DASInstallCompletion(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory(prefix="das_completion_",
            dir=os.environ.get("PRIVA_TEST_TMPDIR", os.environ.get("TMPDIR")))
        self.addCleanup(self.temporary.cleanup)
        self.folder = Path(self.temporary.name)
        self.cache = self.folder / "managed cache"
        self.cache.mkdir()
        self.am = self.folder / "am.vcf.gz"
        self.am.write_bytes(b"opaque external AlphaMissense VCF")
        # A typed statistics stand-in in the fixture priva's format.
        self.stats = self.cache / "input.stats.parquet"
        self.stats.write_text(json.dumps({"command": "alphamissense-domain", "stage": "final",
                                          "vcf": str(self.am),
                                          "builder_source_sha256": builder("alphamissense-domain")}))
        self.mapping = self.folder / "selected mapping.pkl.gz"
        self.mapping.write_bytes(b"opaque fixture mapping")
        self.config_path = self.folder / "config.json"
        self.config = {"assembly": "GRCh38", "threads": 1,
            "alphamissense_pd_stat": str(self.stats),
            "alphamissense_vep_vcf": str(self.am),
            "interpro_mapping_pickle": str(self.mapping)}
        self.config_path.write_text(json.dumps(self.config))
        self.receipt = self.cache / "domain_analysis.hg38.install.json"
        self.priva = fake_priva_binary.environment(self.folder)
        self.history = self.cache / "input.stats.exon_history.seed0.parquet"
        # Built intermediates must be newer than their sources at any timestamp resolution.
        for path in (self.am, self.stats):
            os.utime(path, ns=(10 ** 18, 10 ** 18))


    def native(self):
        return fake_priva_binary.calls(self.priva["PRIVA_FAKE_CALLS"])


    def calls(self):
        """The das-domains analyze calls, oldest first."""
        return [call for call in self.native() if call[:2] == ["das-domains", "analyze"]]


    def test_receipt_device_portability_and_remaining_identity_strictness(self):
        output = self.cache / 'portable-output'
        output.write_bytes(b'output')
        args = SimpleNamespace(profile='hg38', input=[str(self.stats)], output=[str(output)],
                               receipt=str(self.receipt))
        args.expected_inputs = json.dumps(STATE.input_state(args))
        STATE.record(args)
        original = self.receipt.read_bytes()
        real_stat = Path.stat
        def altered(field, target):
            def read(path, *a, **kw):
                value = real_stat(path, *a, **kw)
                attrs = {name: getattr(value, name) for name in dir(value) if name.startswith('st_')}
                if path == target:
                    attrs[field] += 1
                return SimpleNamespace(**attrs)
            return read
        for target in (self.stats, output):
            with self.subTest(target=target.name, field='st_dev'), \
                    patch.object(Path, 'stat', altered('st_dev', target)):
                STATE.check(args)
            for field in ('st_ino', 'st_size', 'st_mtime_ns'):
                with self.subTest(target=target.name, field=field), \
                        patch.object(Path, 'stat', altered(field, target)), self.assertRaises(ValueError):
                    STATE.check(args)
        for group in ('inputs', 'outputs'):
            for field in ('path', 'target'):
                value = json.loads(original)
                value[group][0][field] += '.different'
                self.receipt.write_text(json.dumps(value))
                with self.subTest(group=group, field=field), self.assertRaises(ValueError):
                    STATE.check(args)
        old = json.loads(original)
        old['version'] = 1
        self.receipt.write_text(json.dumps(old))
        with self.assertRaises(ValueError): STATE.check(args)
        self.receipt.write_bytes(original)
        self.assertIn('node_local_st_dev_excluded', json.loads(original)['identity_policy'])
        self.assertNotIn('device', json.loads(original)['inputs'][0])
        STATE.check(args)

    def state_fixture(self):
        source = self.folder / 'analysis.py'; source.write_bytes(b'original source content\n')
        bundle = self.folder / 'bundled.snapshot'; bundle.write_bytes(b'opaque nonempty bundled snapshot')
        output = self.cache / 'portable-output'; output.write_bytes(b'output')
        args = SimpleNamespace(profile='hg38', input=[str(self.stats)], output=[str(output)],
                               source=['analysis=' + str(source)], bundled=['curated=' + str(bundle)],
                               parameter=['filter=enabled', 'threshold=0.6'], receipt=str(self.receipt))
        args.expected_inputs = json.dumps(STATE.input_state(args)); STATE.record(args)
        return args, source, bundle

    def test_source_relocation_with_identical_content_reuses(self):
        args, source, _ = self.state_fixture()
        relocated = self.folder / 'new checkout' / 'renamed.py'; relocated.parent.mkdir()
        relocated.write_bytes(source.read_bytes()); args.source = ['analysis=' + str(relocated)]
        STATE.check(args)
        self.assertNotEqual(source.stat().st_ino, relocated.stat().st_ino)

    def test_changed_source_content_with_same_size_mtime_invalidates(self):
        args, source, _ = self.state_fixture(); before = source.stat()
        source.write_bytes(b'CHANGED! source content\n')
        self.assertEqual(source.stat().st_size, before.st_size)
        os.utime(source, ns=(before.st_atime_ns, before.st_mtime_ns))
        with self.assertRaises(ValueError): STATE.check(args)

    def test_changed_source_before_record_never_promotes_old_receipt(self):
        args, source, _ = self.state_fixture(); old = self.receipt.read_bytes()
        source.write_bytes(b'changed after capture')
        with self.assertRaisesRegex(ValueError, 'changed during generation'): STATE.record(args)
        self.assertEqual(self.receipt.read_bytes(), old)

    def test_source_role_and_missing_source_invalidate(self):
        args, source, _ = self.state_fixture()
        args.source = ['another-role=' + str(source)]
        with self.assertRaises(ValueError): STATE.check(args)
        args.source = ['analysis=' + str(source)]; source.unlink()
        with self.assertRaises(OSError): STATE.check(args)

    def test_identical_relocated_scientific_data_is_still_a_new_selection(self):
        args, _, _ = self.state_fixture(); other = self.folder / 'other-data.pkl'
        other.write_bytes(self.stats.read_bytes()); args.input = [str(other)]
        with self.assertRaises(ValueError): STATE.check(args)

    def test_changed_scientific_input_and_filter_parameter_invalidate(self):
        args, _, _ = self.state_fixture()
        self.stats.write_bytes(b'changed scientific input')
        with self.assertRaises(ValueError): STATE.check(args)
        args.expected_inputs = json.dumps(STATE.input_state(args)); STATE.record(args)
        args.parameter = ['filter=disabled', 'threshold=0.6']
        with self.assertRaises(ValueError): STATE.check(args)
        args.parameter = ['filter=enabled', 'threshold=0.7']
        with self.assertRaises(ValueError): STATE.check(args)
        args.parameter = ['threshold=0.6', 'filter=enabled']
        with self.assertRaises(ValueError): STATE.check(args)

    def test_bundled_resource_uses_only_nonempty_existence(self):
        args, _, bundle = self.state_fixture()
        relocated = self.folder / 'other-bundle'; relocated.write_bytes(b'different nonempty opaque bytes')
        args.bundled = ['curated=' + str(relocated)]
        original_open = Path.open
        def opened(path, *a, **kw):
            if path in (bundle, relocated): raise AssertionError('Bundled content must never be opened')
            return original_open(path, *a, **kw)
        with patch.object(Path, 'open', opened): STATE.check(args)
        relocated.write_bytes(b'')
        with self.assertRaisesRegex(ValueError, 'bundled'): STATE.check(args)
        relocated.unlink()
        with self.assertRaisesRegex(ValueError, 'bundled'): STATE.check(args)

    def test_old_metadata_only_receipt_cannot_be_promoted(self):
        args, _, _ = self.state_fixture(); old = json.loads(self.receipt.read_text())
        old['version'] = 2; self.receipt.write_text(json.dumps(old))
        with self.assertRaises(ValueError): STATE.check(args)
        self.assertEqual(json.loads(self.receipt.read_text())['version'], 2)


if __name__ == "__main__":
    unittest.main()
