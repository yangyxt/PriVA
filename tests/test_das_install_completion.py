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
PRELUDE = r"""
source "$1"
DATA_DIR="$TEST_FOLDER/bundled data"
read_yaml() {
    "$TEST_PYTHON" -c 'import json,sys; print(json.load(open(sys.argv[1])).get(sys.argv[2]) or "null")' "$1" "$2"
}
update_yaml() {
    [[ $2 != "$FAIL_KEY" ]] || return 71
    "$TEST_PYTHON" -c 'import json,sys; p,k,v=sys.argv[1:]; d=json.load(open(p)); d[k]=v; open(p,"w").write(json.dumps(d))' "$@"
}
# The fixture VCF is opaque; its header declares two contigs.
bcftools() { printf '##fileformat=VCFv4.2\n##contig=<ID=chr1,length=10>\n##contig=<ID=chr2,length=10>\n'; }
# The fixture priva runs das_install_state.py through this exported function.
python() { "$TEST_PYTHON" "$@"; }
export -f python
"""


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

    def run_install(self, fail="", fail_key="none", change=None, override="",
                    function="AlphaMissense_pick_intolerant_domains", **extra):
        env = dict(os.environ, TEST_FOLDER=str(self.folder), TEST_PYTHON=sys.executable,
                   FAIL_KEY=fail_key, PYTHONDONTWRITEBYTECODE="1", PRIVA_FAKE_FAIL=fail,
                   # Run the builder as main_install does inside its private view;
                   # native .rkyv values would otherwise start a private installation.
                   PRIVA_CACHE_INSTALL_INTERNAL="1", **self.priva, **extra)
        if change is not None:
            env["PRIVA_FAKE_CHANGE"] = str(change)
        return subprocess.run(["bash", "--noprofile", "--norc", "-c",
            PRELUDE + override + f'\n{function} "$2"',
            "test", str(ROOT / "scripts/install_utils.sh"), str(self.config_path)],
            env=env, text=True, capture_output=True, timeout=60)

    def native(self):
        return fake_priva_binary.calls(self.priva["PRIVA_FAKE_CALLS"])

    def calls(self):
        """The das-domains analyze calls, oldest first."""
        return [call for call in self.native() if call[:2] == ["das-domains", "analyze"]]

    def assert_ok(self, result):
        self.assertEqual(result.returncode, 0, result.stdout + result.stderr)

    def test_cold_then_warm_reuse_and_selected_interpro(self):
        self.assert_ok(self.run_install())
        first = self.receipt.read_bytes()
        calls = self.calls()
        self.assertEqual(len(calls), 1)
        args = calls[0]
        self.assertEqual(args[args.index("--interpro") + 1], str(self.mapping))
        self.assertEqual(args[args.index("--assembly") + 1], "hg38")
        self.assertEqual(self.native()[0], ["cache-install", "resolve", "--key", "interpro_mapping_pickle",
            "--build-config", str(self.config_path), "--source-root", str(ROOT)])
        self.assert_ok(self.run_install(fail="das-domains analyze"))
        self.assertEqual(self.calls(), calls)
        self.assertEqual(self.receipt.read_bytes(), first)
        # The source roles hash the builder identities priva builder-identity prints.
        sources = json.loads(first)["sources"]
        self.assertEqual([entry["role"] for entry in sources], ["analysis", "statistics", "exon-history"])
        identity = json.dumps({"builder": "das-domains", "builder_source_sha256": builder("das-domains")},
                              separators=(",", ":")) + "\n"
        self.assertEqual(sources[0]["sha256"], hashlib.sha256(identity.encode()).hexdigest())

    def test_exon_history_resumes_and_keeps_its_recorded_builder(self):
        shards = self.cache / (self.history.name + ".chromosomes")
        self.assertNotEqual(self.run_install(fail="alpha-exon-history merge").returncode, 0)
        self.assertEqual(sorted(path.name for path in shards.iterdir()), ["00000.parquet", "00001.parquet"])
        self.assertFalse(self.receipt.exists())
        self.assert_ok(self.run_install())
        chromosomes = [call for call in self.native() if call[:2] == ["alpha-exon-history", "chromosome"]]
        self.assertEqual(chromosomes, [
            ["alpha-exon-history", "chromosome", "--vcf", str(self.am), "--chrom", contig,
             "--output", str(shards / name)] for contig, name in (("chr1", "00000.parquet"), ("chr2", "00001.parquet"))])
        self.assertIn(["alpha-exon-history", "merge", "--input", str(shards / "00000.parquet"),
                       "--input", str(shards / "00001.parquet"), "--hash-seed", "0",
                       "--output", str(self.history)], self.native())
        self.assertFalse(shards.exists())
        # A newer priva keeps the recorded history and binds it explicitly.
        count = len(self.native())
        self.assert_ok(self.run_install(PRIVA_FAKE_BUILDER_SHA256="f" * 64))
        later = self.native()[count:]
        self.assertFalse([call for call in later if call[0] == "alpha-exon-history" and
                          call[1] in ("chromosome", "merge", "export")])
        transcripts = next(call for call in later if call[:2] == ["das-domains", "transcripts"])
        self.assertEqual(transcripts[transcripts.index("--history-builder-sha256") + 1], builder("alpha-exon-history"))
        self.assertEqual(len(self.calls()), 2)
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

    def filter_inputs(self):
        spots = self.folder / 'bundled data/HCSeeker/HC_spots.hg38.tsv'
        spots.parent.mkdir(parents=True); spots.write_bytes(b'opaque bundled HCSeeker fixture')
        clinvar = self.folder / 'clinvar.vcf.gz'; clinvar.write_bytes(b'opaque external ClinVar VCF')
        config = json.loads(self.config_path.read_text())
        config.update(clinvar_vcf=str(clinvar))
        self.config_path.write_text(json.dumps(config))
        return spots, self.am, clinvar

    def test_native_argv_and_filter_options(self):
        self.assert_ok(self.run_install())
        self.assertEqual(self.calls()[0], ['das-domains', 'analyze', '--statistics', str(self.stats),
            '--source-builder-sha256', builder('alphamissense-domain'), '--interpro', str(self.mapping),
            '--normalization-json', str(ROOT / 'data/InterPro/InterPro_domain_norm_map.json'),
            '--python-hash-seed', '0', '--threads', '1', '--assembly', 'hg38', '--output-dir', str(self.cache)])
        self.assertFalse([call for call in self.native() if call[0] == 'domain-metadata'])
        self.receipt.unlink()
        spots, am, clinvar = self.filter_inputs(); self.assert_ok(self.run_install())
        metadata = self.cache / 'domain_metadata.hg38.parquet'
        lines = self.cache / 'input.stats.exon_history.seed0.lines'
        native = self.native()
        self.assertIn(['domain-metadata', 'build', '--vcf', str(am), '--output', str(metadata),
                       '--assembly', 'hg38', '--threads', '1'], native)
        self.assertEqual([call for call in native if call[:2] == ['das-domains', 'transcripts']][-1], [
            'das-domains', 'transcripts', '--statistics', str(self.stats), '--exon-history', str(lines),
            '--source-builder-sha256', builder('alphamissense-domain'), '--python-hash-seed', '0',
            '--output', str(self.cache / 'transcript_pp_domain_mapping.hg38.rkyv'),
            '--history-builder-sha256', builder('alpha-exon-history')])
        self.assertEqual(self.calls()[-1], ['das-domains', 'analyze', '--statistics', str(self.stats),
            '--source-builder-sha256', builder('alphamissense-domain'), '--interpro', str(self.mapping),
            '--normalization-json', str(ROOT / 'data/InterPro/InterPro_domain_norm_map.json'),
            '--python-hash-seed', '0', '--threads', '1', '--assembly', 'hg38', '--output-dir', str(self.cache),
            '--metadata', str(metadata), '--spots', str(spots), '--clinvar-vcf', str(clinvar)])
        config = json.loads(self.config_path.read_text())
        self.assertEqual(config['all_intolerant_domains'], str(self.cache / 'all_intolerant_domains.hg38.rkyv'))
        self.assertEqual(config['alphamissense_tranx_domain_map'],
                         str(self.cache / 'transcript_pp_domain_mapping.hg38.rkyv'))

    def test_seed_setting_reaches_every_seeded_builder(self):
        config = json.loads(self.config_path.read_text())
        config['priva_python_hash_seed'] = 7
        self.config_path.write_text(json.dumps(config))
        self.assert_ok(self.run_install())
        native = self.native()
        merge = next(call for call in native if call[:2] == ['alpha-exon-history', 'merge'])
        self.assertEqual(merge[merge.index('--hash-seed') + 1], '7')
        self.assertEqual(merge[-1], str(self.cache / 'input.stats.exon_history.seed7.parquet'))
        for name in ('transcripts', 'analyze'):
            call = next(call for call in native if call[:2] == ['das-domains', name])
            self.assertEqual(call[call.index('--python-hash-seed') + 1], '7')

    def test_filter_enablement_and_external_filter_data_changes_rebuild(self):
        spots, am, clinvar = self.filter_inputs(); self.assert_ok(self.run_install())
        self.assert_ok(self.run_install(fail='das-domains analyze')); self.assertEqual(len(self.calls()), 1)
        am.write_bytes(b'changed selected AM VCF'); self.assert_ok(self.run_install())
        self.assertEqual(len(self.calls()), 2)
        clinvar.write_bytes(b'changed selected ClinVar VCF'); self.assert_ok(self.run_install())
        self.assertEqual(len(self.calls()), 3)
        spots.unlink(); self.assert_ok(self.run_install()); self.assertEqual(len(self.calls()), 4)
        clinvar.unlink(); self.assert_ok(self.run_install()); self.assertEqual(len(self.calls()), 5)

    def test_selected_interpro_archive_is_a_strict_input(self):
        self.assert_ok(self.run_install())
        self.mapping.write_bytes(b'changed selected mapping'); self.assert_ok(self.run_install())
        self.assertEqual(len(self.calls()), 2)

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
            self.assert_ok(self.run_install(fail='das-domains analyze'))
            self.assertEqual(json.loads(self.config_path.read_text()), ready)
            self.assertEqual(len(self.calls()), count)
        # Already normalized values must not call update_yaml at all.
        for key in keys:
            config_before = (self.config_path.read_bytes(), self.config_path.stat().st_mtime_ns)
            self.assert_ok(self.run_install(fail='das-domains analyze', fail_key=key))
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
            self.assertNotEqual(self.run_install(fail='das-domains analyze', fail_key=failure).returncode, 0)
            self.assertEqual(self.receipt.read_bytes(), receipt)
            self.assertEqual(len(self.calls()), count)
            self.assert_ok(self.run_install(fail='das-domains analyze'))
            self.assertEqual(json.loads(self.config_path.read_text()), ready)
            self.assertEqual(len(self.calls()), count)

    def test_config_path_with_spaces(self):
        selected = self.folder / "config with spaces.json"
        self.config_path.rename(selected)
        self.config_path = selected
        self.assert_ok(self.run_install())
        self.assertTrue(self.receipt.is_file())

    def test_failed_native_stages_never_reuse_partial_generation(self):
        # The fixture records a failing call before it exits.
        for stage, failed_analyses in (("das-domains transcripts", 0), ("das-domains analyze", 1),
                                       ("das-install-state record", 1)):
            with self.subTest(stage=stage):
                self.assert_ok(self.run_install())
                self.mapping.write_bytes(self.mapping.read_bytes() + b"!")
                count = len(self.calls())
                self.assertNotEqual(self.run_install(fail=stage).returncode, 0)
                self.assertFalse(self.receipt.exists())
                self.assert_ok(self.run_install())
                self.assertEqual(len(self.calls()), count + failed_analyses + 1)
                self.assertTrue(self.receipt.is_file())

    def test_yaml_failures_and_selected_destination_changes_retry(self):
        for key in ("alphamissense_intolerant_domains", "alphamissense_tranx_domain_map", "all_intolerant_domains"):
            with self.subTest(key=key):
                config = json.loads(self.config_path.read_text())
                suffix = ".tsv" if key == "alphamissense_intolerant_domains" else ".rkyv"
                config[key] = str(self.cache / (key + ".selected" + suffix))
                self.config_path.write_text(json.dumps(config))
                self.assertNotEqual(self.run_install(fail_key=key).returncode, 0)
                self.assertFalse(self.receipt.exists())
                count = len(self.calls())
                self.assert_ok(self.run_install())
                self.assertEqual(len(self.calls()), count + 1)
                self.assertTrue(Path(config[key]).is_file())

    def test_selected_pickle_destinations_select_native_siblings(self):
        config = json.loads(self.config_path.read_text())
        config.update(alphamissense_tranx_domain_map=str(self.cache / "selected map.pkl"),
                      all_intolerant_domains=str(self.cache / "selected domains.pkl.gz"))
        self.config_path.write_text(json.dumps(config))
        self.assert_ok(self.run_install())
        config = json.loads(self.config_path.read_text())
        self.assertEqual(config["alphamissense_tranx_domain_map"], str(self.cache / "selected map.rkyv"))
        self.assertEqual(config["all_intolerant_domains"], str(self.cache / "selected domains.rkyv"))
        self.assertTrue(Path(config["all_intolerant_domains"]).is_file())

    def test_changed_input_during_command_cannot_record_completion(self):
        self.assertNotEqual(self.run_install(change=self.mapping).returncode, 0)
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

    def legacy_pickle(self):
        stats = self.cache / "input.stats.pkl"
        stats.write_bytes(b"completed Python statistics fixture")
        config = json.loads(self.config_path.read_text())
        config["alphamissense_pd_stat"] = str(stats)
        self.config_path.write_text(json.dumps(config))
        return stats

    def test_legacy_python_receipt_keeps_python_outputs_without_native_calls(self):
        self.legacy_pickle()
        # Stand in for a valid Python completion receipt of an existing installation.
        override = r"""
python() {
    if [[ $1 == */das_install_state.py && $2 == check ]]; then
        printf 'LEGACY_CHECK %s\n' "${@:3}" >> "$TEST_FOLDER/legacy.events"; return 0
    fi
    "$TEST_PYTHON" "$@"
}
"""
        for function in ("AlphaMissense_stat", "AlphaMissense_pick_intolerant_domains"):
            with self.subTest(function=function):
                self.assert_ok(self.run_install(override=override, function=function))
        self.assertEqual(self.native(), [])
        events = (self.folder / "legacy.events").read_text()
        self.assertIn("analysis=" + str(ROOT / "scripts/am_pick_intolerant_domains.py"), events)
        config = json.loads(self.config_path.read_text())
        self.assertEqual(config["alphamissense_pd_stat"], str(self.cache / "input.stats.pkl"))
        self.assertEqual(config["all_intolerant_domains"], str(self.cache / "all_intolerant_domains.hg38.pkl"))
        self.assertEqual(config["alphamissense_tranx_domain_map"],
                         str(self.cache / "transcript_pp_domain_mapping.hg38.pkl"))

    def test_legacy_pickle_without_python_receipt_builds_native_statistics(self):
        self.legacy_pickle()
        result = self.run_install()
        self.assertNotEqual(result.returncode, 0)
        self.assertIn("run AlphaMissense_stat", result.stderr)
        self.assertEqual(self.native(), [])
        self.stats.unlink()
        self.assert_ok(self.run_install(function="AlphaMissense_stat"))
        shards = self.cache / "input.stats.parquet.chromosomes"
        self.assertEqual([call for call in self.native() if call[0] == "alphamissense-domain"], [
            ["alphamissense-domain", "chromosome", "--vcf", str(self.am), "--chrom", "chr1",
             "--output", str(shards / "00000.parquet")],
            ["alphamissense-domain", "chromosome", "--vcf", str(self.am), "--chrom", "chr2",
             "--output", str(shards / "00001.parquet")],
            ["alphamissense-domain", "merge", "--input", str(shards / "00000.parquet"),
             "--input", str(shards / "00001.parquet"), "--output", str(self.stats)],
            ["alphamissense-domain", "check", "--input", str(self.stats), "--vcf", str(self.am)]])
        self.assertEqual(json.loads(self.config_path.read_text())["alphamissense_pd_stat"], str(self.stats))
        # A complete Parquet keeps its recorded builder.
        count = len(self.native())
        self.assert_ok(self.run_install(function="AlphaMissense_stat", PRIVA_FAKE_BUILDER_SHA256="f" * 64))
        self.assertEqual(self.native()[count:], [
            ["alphamissense-domain", "identity", "--input", str(self.stats)],
            ["alphamissense-domain", "check", "--input", str(self.stats), "--vcf", str(self.am),
             "--builder-sha256", builder("alphamissense-domain")]])
        self.assert_ok(self.run_install())
        self.assertEqual(len(self.calls()), 1)


if __name__ == "__main__":
    unittest.main()
