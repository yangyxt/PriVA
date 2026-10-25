"""Real installer/state helper, synthetic scientific writer and injected failures.

These tests establish orchestration/recovery, not domain model acceptance.
"""
import json
import os
from pathlib import Path
import re
import subprocess
import sys
import tempfile
import unittest
from types import SimpleNamespace
from unittest.mock import patch
import importlib.util


ROOT = Path(__file__).resolve().parents[1]
_spec = importlib.util.spec_from_file_location('das_install_state_test', ROOT / 'scripts/das_install_state.py')
STATE = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(STATE)
PRELUDE = r'''
source "$1"
DATA_DIR="$TEST_FOLDER/bundled data"
read_yaml() {
    "$TEST_PYTHON" -c 'import json,sys; print(json.load(open(sys.argv[1])).get(sys.argv[2]) or "null")' "$1" "$2"
}
update_yaml() {
    [[ $2 != "$FAIL_KEY" ]] || return 71
    "$TEST_PYTHON" -c 'import json,sys; p,k,v=sys.argv[1:]; d=json.load(open(p)); d[k]=v; open(p,"w").write(json.dumps(d))' "$@"
}
python() {
    if [[ $1 == */am_pick_intolerant_domains.py ]]; then
        "$TEST_PYTHON" "$TEST_FOLDER/writer.py" "$@"
    else
        "$TEST_PYTHON" "$@"
    fi
}
'''
WRITER = r'''
import json,os,sys
from pathlib import Path
args=sys.argv[1:]
with open(os.environ["TEST_FOLDER"]+"/calls.jsonl","a") as f:
    f.write(json.dumps(args)+"\n")
folder=Path(args[args.index("--output_dir")+1])
mode=os.environ["FAIL_MODE"]
for index,(stem,suffix) in enumerate((("transcript_pp_domain_mapping","pkl"),
        ("domain_tolerance_analysis","tsv"),("all_intolerant_domains","pkl"))):
    path=folder/f"{stem}.hg38.{suffix}"
    path.write_bytes(b"complete synthetic output:"+stem.encode())
    if mode == str(index):
        raise SystemExit(73)
if mode == "empty": path.write_bytes(b"")
if mode == "input_changed":
    Path(args[args.index("--interpro_mapping_pickle")+1]).write_bytes(b"changed selected input")
'''


class DASInstallCompletion(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory(prefix="das_completion_",
            dir=os.environ.get("PRIVA_TEST_TMPDIR", os.environ.get("TMPDIR")))
        self.addCleanup(self.temporary.cleanup)
        self.folder = Path(self.temporary.name)
        self.cache = self.folder / "managed cache"
        self.cache.mkdir()
        self.stats = self.cache / "input.stats.pkl"
        self.stats.write_bytes(b"opaque fixture input")
        self.mapping = self.folder / "selected mapping.pkl.gz"
        self.mapping.write_bytes(b"opaque fixture mapping")
        self.config_path = self.folder / "config.json"
        self.config = {"assembly": "GRCh38", "threads": 1,
            "alphamissense_pd_stat": str(self.stats),
            "interpro_mapping_pickle": str(self.mapping)}
        self.config_path.write_text(json.dumps(self.config))
        (self.folder / "writer.py").write_text(WRITER)
        self.receipt = self.cache / "domain_analysis.hg38.install.json"

    def run_install(self, mode="none", fail_key="none", override="", installer=None):
        env = dict(os.environ, TEST_FOLDER=str(self.folder), TEST_PYTHON=sys.executable,
                   FAIL_MODE=mode, FAIL_KEY=fail_key, PYTHONDONTWRITEBYTECODE="1")
        return subprocess.run(["bash", "--noprofile", "--norc", "-c",
            PRELUDE + override + '\nAlphaMissense_pick_intolerant_domains "$2"',
            "test", str(installer or ROOT / "scripts/install_utils.sh"), str(self.config_path)],
            env=env, text=True, capture_output=True, timeout=60)

    def calls(self):
        path = self.folder / "calls.jsonl"
        return [json.loads(line) for line in path.read_text().splitlines()] if path.exists() else []

    def assert_ok(self, result):
        self.assertEqual(result.returncode, 0, result.stdout + result.stderr)

    def test_cold_then_warm_reuse_and_selected_interpro(self):
        self.assert_ok(self.run_install())
        first = self.receipt.read_bytes()
        calls = self.calls()
        self.assertEqual(len(calls), 1)
        args = calls[0]
        self.assertEqual(args[args.index("--interpro_mapping_pickle") + 1], str(self.mapping))
        self.assertEqual(args[args.index("--assembly") + 1], "hg38")
        self.assert_ok(self.run_install(mode="0"))
        self.assertEqual(self.calls(), calls)
        self.assertEqual(self.receipt.read_bytes(), first)

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

    def test_installer_reuses_identical_sources_from_other_checkout(self):
        peer = ROOT.parent / ('PriVA-nextflow' if ROOT.name == 'PriVA' else 'PriVA')
        installer = peer / 'scripts/install_utils.sh'
        if not installer.is_file(): self.skipTest('A second installer checkout is required')
        self.assert_ok(self.run_install()); before = self.receipt.read_bytes()
        self.assert_ok(self.run_install(mode='0', installer=installer))
        self.assertEqual(len(self.calls()), 1)
        self.assertEqual(self.receipt.read_bytes(), before)

    def filter_inputs(self):
        spots = self.folder / 'bundled data/HCSeeker/HC_spots.hg38.tsv'
        spots.parent.mkdir(parents=True); spots.write_bytes(b'opaque bundled HCSeeker fixture')
        am = self.folder / 'am.vcf.gz'; am.write_bytes(b'opaque external AlphaMissense VCF')
        clinvar = self.folder / 'clinvar.vcf.gz'; clinvar.write_bytes(b'opaque external ClinVar VCF')
        config = json.loads(self.config_path.read_text())
        config.update(alphamissense_vep_vcf=str(am), clinvar_vcf=str(clinvar))
        self.config_path.write_text(json.dumps(config))
        return spots, am, clinvar

    def test_scientific_argv_order_and_options_remain_unchanged(self):
        spots, am, clinvar = self.filter_inputs(); self.assert_ok(self.run_install())
        self.assertEqual(self.calls()[0], [str(ROOT / 'scripts/am_pick_intolerant_domains.py'),
            '--assembly', 'hg38', '--pickle_file', str(self.stats), '--interpro_mapping_pickle', str(self.mapping),
            '--threads', '1', '--output_dir', str(self.cache), '--am_vcf', str(am),
            '--hcseeker_spots_tsv', str(spots), '--clinvar_vcf', str(clinvar)])

    def test_filter_enablement_and_external_filter_data_changes_rebuild(self):
        spots, am, clinvar = self.filter_inputs(); self.assert_ok(self.run_install())
        self.assert_ok(self.run_install(mode='0')); self.assertEqual(len(self.calls()), 1)
        am.write_bytes(b'changed selected AM VCF'); self.assert_ok(self.run_install())
        self.assertEqual(len(self.calls()), 2)
        clinvar.write_bytes(b'changed selected ClinVar VCF'); self.assert_ok(self.run_install())
        self.assertEqual(len(self.calls()), 3)
        spots.unlink(); self.assert_ok(self.run_install()); self.assertEqual(len(self.calls()), 4)
        clinvar.unlink(); self.assert_ok(self.run_install()); self.assertEqual(len(self.calls()), 5)

    def test_only_actual_bundled_mapping_gets_existence_only_reuse(self):
        bundled = self.folder / 'bundled data/InterPro/Interpro_entry_mapping.pkl.gz'
        bundled.parent.mkdir(parents=True); bundled.write_bytes(b'opaque bundled mapping')
        self.config['interpro_mapping_pickle'] = str(bundled)
        self.config_path.write_text(json.dumps(self.config)); self.assert_ok(self.run_install())
        bundled.write_bytes(b'changed but nonempty bundled mapping')
        self.assert_ok(self.run_install(mode='0')); self.assertEqual(len(self.calls()), 1)
        external = self.folder / 'Interpro_entry_mapping.pkl.gz'; external.write_bytes(b'external selected mapping')
        config = json.loads(self.config_path.read_text()); config['interpro_mapping_pickle'] = str(external)
        self.config_path.write_text(json.dumps(config)); self.assert_ok(self.run_install())
        external.write_bytes(b'changed external mapping'); self.assert_ok(self.run_install())
        self.assertEqual(len(self.calls()), 3)

    def test_bundled_mapping_link_to_external_remains_strict(self):
        bundled = self.folder / 'bundled data/InterPro/Interpro_entry_mapping.pkl.gz'
        bundled.parent.mkdir(parents=True); bundled.symlink_to(self.mapping)
        self.config['interpro_mapping_pickle'] = str(bundled)
        self.config_path.write_text(json.dumps(self.config)); self.assert_ok(self.run_install())
        self.mapping.write_bytes(b'changed external symlink target'); self.assert_ok(self.run_install())
        self.assertEqual(len(self.calls()), 2)

    def test_warm_default_keys_repaired_without_scientific_or_cache_write(self):
        self.assert_ok(self.run_install())
        ready = json.loads(self.config_path.read_text())
        keys = ('alphamissense_intolerant_domains', 'alphamissense_tranx_domain_map', 'all_intolerant_domains')
        selected = [Path(ready[key]) for key in keys] + [self.receipt]
        before = {p: (p.read_bytes(), p.stat().st_ino, p.stat().st_mtime_ns) for p in selected}
        count = len(self.calls())
        for missing in ('absent', 'blank', 'null'):
            value = dict(ready)
            for key in keys:
                if missing == 'absent': value.pop(key)
                else: value[key] = '' if missing == 'blank' else None
            self.config_path.write_text(json.dumps(value))
            self.assert_ok(self.run_install(mode='0'))
            self.assertEqual(json.loads(self.config_path.read_text()), ready)
            self.assertEqual(len(self.calls()), count)
        # Already normalized values must not call update_yaml at all.
        for key in keys:
            config_before = (self.config_path.read_bytes(), self.config_path.stat().st_mtime_ns)
            self.assert_ok(self.run_install(mode='0', fail_key=key))
            self.assertEqual(config_before, (self.config_path.read_bytes(), self.config_path.stat().st_mtime_ns))
        self.assertEqual(before, {p: (p.read_bytes(), p.stat().st_ino, p.stat().st_mtime_ns) for p in selected})

    def test_warm_default_update_failure_propagates_and_retry_reuses(self):
        self.assert_ok(self.run_install())
        ready = json.loads(self.config_path.read_text())
        keys = ('alphamissense_intolerant_domains', 'alphamissense_tranx_domain_map', 'all_intolerant_domains')
        count = len(self.calls())
        receipt = self.receipt.read_bytes()
        for failure in keys:
            value = dict(ready, **{key: None for key in keys})
            self.config_path.write_text(json.dumps(value))
            self.assertNotEqual(self.run_install(mode='0', fail_key=failure).returncode, 0)
            self.assertEqual(self.receipt.read_bytes(), receipt)
            self.assertEqual(len(self.calls()), count)
            self.assert_ok(self.run_install(mode='0'))
            self.assertEqual(json.loads(self.config_path.read_text()), ready)
            self.assertEqual(len(self.calls()), count)

    def test_config_path_with_spaces(self):
        selected = self.folder / "config with spaces.json"
        self.config_path.rename(selected)
        self.config_path = selected
        self.assert_ok(self.run_install())
        self.assertTrue(self.receipt.is_file())

    def test_failed_scientific_stages_never_reuse_partial_generation(self):
        for stage in ("0", "1", "2", "empty"):
            with self.subTest(stage=stage):
                self.assert_ok(self.run_install())
                self.mapping.write_bytes(self.mapping.read_bytes() + b"!")
                count = len(self.calls())
                self.assertNotEqual(self.run_install(mode=stage).returncode, 0)
                self.assertFalse(self.receipt.exists())
                self.assert_ok(self.run_install())
                self.assertEqual(len(self.calls()), count + 2)
                self.assertTrue(self.receipt.is_file())

    def test_yaml_failures_and_selected_destination_changes_retry(self):
        for key in ("alphamissense_intolerant_domains", "alphamissense_tranx_domain_map", "all_intolerant_domains"):
            with self.subTest(key=key):
                config = json.loads(self.config_path.read_text())
                config[key] = str(self.cache / (key + ".selected"))
                self.config_path.write_text(json.dumps(config))
                self.assertNotEqual(self.run_install(fail_key=key).returncode, 0)
                self.assertFalse(self.receipt.exists())
                count = len(self.calls())
                self.assert_ok(self.run_install())
                self.assertEqual(len(self.calls()), count + 1)
                self.assertTrue(Path(config[key]).is_file())

    def test_changed_input_during_command_cannot_record_completion(self):
        self.assertNotEqual(self.run_install(mode="input_changed").returncode, 0)
        self.assertFalse(self.receipt.exists())
        self.assert_ok(self.run_install())

    def test_output_or_receipt_damage_rebuilds(self):
        self.assert_ok(self.run_install())
        config = json.loads(self.config_path.read_text())
        Path(config["all_intolerant_domains"]).write_bytes(b"truncated")
        self.assert_ok(self.run_install())
        self.assertEqual(len(self.calls()), 2)
        self.receipt.write_text("invalid JSON")
        self.assert_ok(self.run_install())
        self.assertEqual(len(self.calls()), 3)

    def test_old_false_reuse_reproduced_and_current_requires_completion(self):
        previous = os.environ.get("PRIVA_TEST_BEFORE_SOURCE_ROOT")
        if not previous:
            self.skipTest("Frozen predecessor required for historical defect reproduction")
        source = (Path(previous) / "scripts/install_utils.sh").read_text()
        old = re.search(r"^function AlphaMissense_pick_intolerant_domains\(\).*?(?=^function )",
                        source, re.M | re.S).group()
        # All three files can look current even when the real command failed
        # immediately after its last main output and before later outputs.
        for stem, suffix in (("domain_tolerance_analysis", "tsv"),
                             ("transcript_pp_domain_mapping", "pkl"),
                             ("all_intolerant_domains", "pkl")):
            path = self.cache / f"{stem}.hg38.{suffix}"
            path.write_bytes(b"partial failed generation")
            future = max(self.stats.stat().st_mtime_ns,
                (ROOT / "scripts/am_pick_intolerant_domains.py").stat().st_mtime_ns) + 1000000000
            os.utime(path, ns=(future, future))
        self.assert_ok(self.run_install(mode="0", override=old))
        self.assertEqual(self.calls(), [])
        self.assert_ok(self.run_install())
        self.assertEqual(len(self.calls()), 1)
        self.assertTrue(self.receipt.is_file())


if __name__ == "__main__":
    unittest.main()
