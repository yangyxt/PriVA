#!/usr/bin/env python3
"""Retained tiny publication fixtures; no real pickle load or calculation.

Use --fixture-root in an explicitly authorized runtime/data tree. Tests retain
their inputs, interrupted files and sidecars. Native execution records below
are labelled synthetic; they provide no evidence about a real AlphaMissense
run. Both installers run this same mirrored suite against their own shell.
"""
from __future__ import annotations

import argparse
import copy
from datetime import datetime, timezone
import hashlib
import importlib.util
import json
import os
from pathlib import Path
import pickle
import subprocess
import sys
import unittest
from unittest.mock import patch

sys.dont_write_bytecode = True
ROOT = Path(__file__).absolute().parents[1]
HELPER = ROOT / 'scripts/alphamissense_original_cache.py'
spec = importlib.util.spec_from_file_location('original_cache_under_test', HELPER)
m = importlib.util.module_from_spec(spec); spec.loader.exec_module(m)
FIXTURE_ROOT = None


class OriginalCache(unittest.TestCase):
    def setUp(self):
        self.folder = FIXTURE_ROOT / self._testMethodName
        self.folder.mkdir()
        self.cache = self.folder / 'cache root'; self.cache.mkdir()
        self.vcf = self.cache / 'input.vcf.gz'; self.vcf.write_bytes(b'SYNTHETIC indexed-input metadata fixture')
        self.index = Path(str(self.vcf) + '.tbi'); self.index.write_bytes(b'SYNTHETIC index')
        self.original = self.folder / 'completed original.pkl'
        self.original.write_bytes(pickle.dumps({'synthetic': [1, 2, 3]}, protocol=5))
        self.source = self.folder / 'producer.py'
        self.source.write_bytes((ROOT / 'scripts/stat_protein_domain_amscores.py').read_bytes())
        self.output = self.cache / 'selected.pkl'
        self.run_record = self.folder / 'synthetic-native-record.txt'
        self.run_record.write_text('SYNTHETIC fixture: COMPLETED 0:0; no scientific execution\n')
        self.run = self.folder / 'synthetic-run.json'
        self.attestation = dict(schema=m.RUN_SCHEMA, complete=True,
            recipe=m.snapshot(source_script=self.source, input_vcf=self.vcf),
            source_binding=dict(basis='recorded-at-run', description='SYNTHETIC pre-run fixture hash'),
            input_binding=dict(basis='retrospective-attestation', description='SYNTHETIC retained input metadata',
                               limitation='Fixture only; no actual producer execution'),
            native_execution=dict(state='COMPLETED', exit_code='0:0',
                                  records=[dict(path=str(self.run_record), sha256=m.digest(self.run_record))]),
            output=dict(path=str(self.original), size=str(self.original.stat().st_size), sha256=m.digest(self.original)))
        self.write_attestation()
        self.kwargs = dict(source_script=self.source, input_vcf=self.vcf, original_pickle=self.original,
                           execution_evidence=self.run, output_pickle=self.output, cache_root=self.cache)

    def write_attestation(self):
        self.run.write_text(json.dumps(self.attestation))

    @property
    def ready(self):
        return Path(str(self.output) + m.READY_SUFFIX)

    @property
    def pending(self):
        return Path(str(self.output) + m.PENDING_SUFFIX)

    def imported(self):
        return m.import_original(**self.kwargs)

    def check(self, **updates):
        options = {key: self.kwargs[key] for key in ('source_script', 'input_vcf', 'output_pickle', 'cache_root')}
        options.update(updates)
        return m.check(**options)

    def interrupted(self, point):
        original_json = m.write_json
        if point == 'receipt':
            def writer(path, value, **kwargs):
                if Path(path) == self.ready:
                    raise OSError('SYNTHETIC interrupted before receipt')
                return original_json(path, value, **kwargs)
            return patch.object(m, 'write_json', writer)
        original_link = m.os.link
        def link(source, target, **kwargs):
            if Path(target) == self.output:
                raise OSError('SYNTHETIC interrupted before payload publication')
            return original_link(source, target, **kwargs)
        return patch.object(m.os, 'link', link)

    def shell(self, *, output=None):
        config = self.folder / 'config.json'
        config.write_text(json.dumps(dict(cache_root=str(self.cache), alphamissense_pd_stat=str(output or self.output),
                                          alphamissense_vep_vcf=str(self.vcf))))
        events = self.folder / 'shell.events'
        script = r'''
source "$TEST_INSTALLER"
log() { printf '%s\n' "$*" >&2; }
read_yaml() { "$TEST_PYTHON" -c 'import json,sys; print(json.load(open(sys.argv[1])).get(sys.argv[2],""))' "$1" "$2"; }
update_yaml() { printf 'UPDATE:%s\n' "$*" >> "$TEST_EVENTS"; }
python() {
    if [[ $1 == */alphamissense_domain_cache.py ]]; then
        printf 'BUILDER:%s\n' "$*" >> "$TEST_EVENTS"
        return 0
    fi
    "$TEST_PYTHON" "$@"
}
AlphaMissense_stat "$TEST_CONFIG"
'''
        env = dict(os.environ, TEST_INSTALLER=str(ROOT / 'scripts/install_utils.sh'),
                   TEST_PYTHON=sys.executable, TEST_CONFIG=str(config), TEST_EVENTS=str(events),
                   PYTHONDONTWRITEBYTECODE='1')
        result = subprocess.run(['bash', '-c', script], env=env, text=True, capture_output=True)
        return result, events.read_text() if events.exists() else ''

    def test_import_preserves_original_and_uses_distinct_receipt(self):
        before = m.entry(self.original, hashed=True)
        result = self.imported()
        self.assertFalse(result['reused'])
        self.assertEqual(self.output.read_bytes(), self.original.read_bytes())
        self.assertEqual(m.entry(self.original, hashed=True), before)
        self.assertFalse(Path(str(self.output) + '.complete.json').exists())
        self.assertFalse(self.pending.exists())
        value = m.read_json(self.ready)
        self.assertEqual(value['schema'], m.SCHEMA)
        self.assertFalse(value['installation_accepted'])
        self.assertEqual(value['execution']['attestation']['input_binding']['limitation'],
                         'Fixture only; no actual producer execution')
        self.assertTrue(self.check()['reused'])

    def test_repeat_import_does_not_change_published_files(self):
        self.imported(); before = [m.entry(path, hashed=True) for path in (self.output, self.ready)]
        self.assertTrue(self.imported()['reused'])
        self.assertEqual([m.entry(path, hashed=True) for path in (self.output, self.ready)], before)

    def test_same_source_content_in_another_checkout_reuses(self):
        self.imported(); other = self.folder / 'other-producer.py'
        other.write_bytes(self.source.read_bytes())
        self.assertTrue(self.check(source_script=other)['reused'])

    def test_warm_check_does_not_hash_the_large_output(self):
        self.imported(); original_digest = m.digest
        def observed(path):
            self.assertNotEqual(Path(path), self.output)
            return original_digest(path)
        with patch.object(m, 'digest', observed):
            self.assertTrue(self.check()['reused'])

    def test_no_original_receipt_returns_unowned(self):
        self.assertIsNone(self.check())
        self.assertFalse(self.output.exists())

    def test_changed_input_blocks_reuse(self):
        self.imported(); self.vcf.write_bytes(b'changed input')
        with self.assertRaisesRegex(ValueError, 'differs'): self.check()

    def test_missing_input_blocks_reuse(self):
        self.imported(); self.vcf.unlink()
        with self.assertRaises(OSError): self.check()

    def test_changed_index_blocks_reuse(self):
        self.imported(); self.index.write_bytes(b'changed index')
        with self.assertRaisesRegex(ValueError, 'differs'): self.check()

    def test_missing_index_blocks_reuse(self):
        self.imported(); self.index.unlink()
        with self.assertRaisesRegex(ValueError, 'indexed'): self.check()

    def test_new_index_changes_selected_recipe(self):
        self.imported(); Path(str(self.vcf) + '.csi').write_bytes(b'new index')
        with self.assertRaisesRegex(ValueError, 'differs'): self.check()

    def test_changed_source_blocks_reuse(self):
        self.imported(); self.source.write_bytes(self.source.read_bytes() + b'\n# changed\n')
        with self.assertRaisesRegex(ValueError, 'pinned'): self.check()

    def test_missing_source_blocks_reuse(self):
        self.imported(); self.source.unlink()
        with self.assertRaises(OSError): self.check()

    def test_changed_output_blocks_reuse_and_import_preserves_it(self):
        self.imported(); self.output.write_bytes(b'changed publication')
        with self.assertRaisesRegex(ValueError, 'changed'): self.check()
        with self.assertRaises(ValueError): self.imported()
        self.assertEqual(self.output.read_bytes(), b'changed publication')

    def test_missing_output_blocks_reuse(self):
        self.imported(); self.output.unlink()
        with self.assertRaises(OSError): self.check()

    def test_receipt_requires_complete_original_native_execution(self):
        self.imported(); value = m.read_json(self.ready)
        value['execution']['attestation']['native_execution']['exit_code'] = '1:0'
        self.ready.write_text(json.dumps(value))
        with self.assertRaisesRegex(ValueError, 'successfully'): self.check()

    def test_builder_receipt_is_not_original_proof(self):
        Path(str(self.output) + '.complete.json').write_text(json.dumps({'schema':'priva-alphamissense-domain-checkpoints-v1'}))
        self.assertIsNone(self.check())
        with self.assertRaisesRegex(ValueError, 'chromosome-builder'): self.imported()
        self.assertFalse(self.output.exists())

    def test_unowned_output_and_metadata_are_preserved(self):
        self.output.write_bytes(b'nonowned output')
        with self.assertRaisesRegex(ValueError, 'Preserve existing output'): self.imported()
        self.assertEqual(self.output.read_bytes(), b'nonowned output')
        self.output.unlink(); self.ready.write_bytes(b'nonowned metadata')
        with self.assertRaises(ValueError): self.imported()
        self.assertEqual(self.ready.read_bytes(), b'nonowned metadata')

    def test_incomplete_or_failed_execution_is_not_adopted(self):
        self.attestation['complete'] = False; self.write_attestation()
        with self.assertRaisesRegex(ValueError, 'completion'): self.imported()
        self.attestation['complete'] = True
        self.attestation['native_execution']['state'] = 'FAILED'; self.write_attestation()
        with self.assertRaisesRegex(ValueError, 'successfully'): self.imported()
        self.assertFalse(self.output.exists())

    def test_retrospective_attestation_must_state_proof_gap(self):
        del self.attestation['input_binding']['limitation']; self.write_attestation()
        with self.assertRaisesRegex(ValueError, 'missing historical proof'): self.imported()

    def test_changed_native_evidence_is_not_adopted(self):
        self.run_record.write_bytes(b'changed native record')
        with self.assertRaisesRegex(ValueError, 'checksum'): self.imported()

    def test_changed_native_evidence_blocks_warm_reuse(self):
        self.imported(); self.run_record.write_bytes(b'changed native record after publication')
        with self.assertRaisesRegex(ValueError, 'execution evidence changed'): self.check()

    def test_changed_run_evidence_blocks_warm_reuse(self):
        self.imported(); self.attestation['input_binding']['description'] += ' changed'
        self.write_attestation()
        with self.assertRaisesRegex(ValueError, 'execution evidence changed'): self.check()

    def test_embedded_attestation_cannot_disagree_with_evidence(self):
        self.imported(); value = m.read_json(self.ready)
        value['execution']['attestation']['input_binding']['description'] += ' forged'
        self.ready.write_text(json.dumps(value))
        with self.assertRaisesRegex(ValueError, 'Embedded'): self.check()

    def test_truncated_original_is_not_published(self):
        self.original.write_bytes(self.original.read_bytes()[:-2])
        with self.assertRaisesRegex(ValueError, 'output attestation'): self.imported()
        self.assertFalse(self.output.exists())

    def test_empty_original_is_not_published(self):
        self.original.write_bytes(b'')
        with self.assertRaisesRegex(ValueError, 'nonempty'): self.imported()

    def test_short_copy_is_not_published_and_retry_retains_it(self):
        def short_copy(original, stage, expected):
            stage.write_bytes(b'short retained bytes')
            raise OSError('SYNTHETIC short write')
        with patch.object(m, 'copy_original', short_copy):
            with self.assertRaisesRegex(OSError, 'short write'): self.imported()
        retained = Path(m.read_json(self.pending)['stage'])
        self.assertEqual(retained.read_bytes(), b'short retained bytes')
        self.assertFalse(self.output.exists())
        self.imported()
        self.assertEqual(retained.read_bytes(), b'short retained bytes')
        self.assertTrue(self.check()['reused'])

    def test_interrupted_before_payload_recovers(self):
        with self.interrupted('payload'):
            with self.assertRaises(OSError): self.imported()
        with self.assertRaisesRegex(ValueError, 'pending'): self.check()
        self.assertFalse(self.output.exists())
        self.imported(); self.assertTrue(self.check()['reused'])

    def test_interrupted_after_payload_recovers_without_recopy(self):
        with self.interrupted('receipt'):
            with self.assertRaises(OSError): self.imported()
        self.assertTrue(self.output.exists())
        with self.assertRaisesRegex(ValueError, 'pending'): self.check()
        with patch.object(m, 'copy_original', side_effect=AssertionError('must not recopy')):
            self.imported()
        self.assertTrue(self.check()['reused'])

    def test_foreign_output_during_recovery_is_preserved(self):
        with self.interrupted('payload'):
            with self.assertRaises(OSError): self.imported()
        self.output.write_bytes(self.original.read_bytes())
        with self.assertRaisesRegex(ValueError, 'not owned'): self.imported()
        self.assertTrue(self.pending.exists())

    def test_changed_input_during_copy_never_publishes(self):
        original_copy = m.copy_original
        def copied(*args):
            result = original_copy(*args)
            self.index.write_bytes(b'changed during publication')
            return result
        with patch.object(m, 'copy_original', copied):
            with self.assertRaisesRegex(ValueError, 'changed during'): self.imported()
        self.assertFalse(self.output.exists()); self.assertFalse(self.ready.exists())

    def test_changed_source_during_copy_never_publishes(self):
        original_copy = m.copy_original
        def copied(*args):
            result = original_copy(*args)
            self.source.write_bytes(b'changed original source')
            return result
        with patch.object(m, 'copy_original', copied):
            with self.assertRaisesRegex(ValueError, 'pinned'): self.imported()
        self.assertFalse(self.output.exists()); self.assertFalse(self.ready.exists())

    def test_input_output_alias_is_rejected(self):
        self.kwargs['output_pickle'] = self.vcf
        with self.assertRaisesRegex(ValueError, 'aliases'): self.imported()

    def test_original_output_hardlink_alias_is_rejected(self):
        os.link(self.original, self.output)
        with self.assertRaisesRegex(ValueError, 'aliases'): self.imported()

    def test_linked_destination_and_parent_escape_are_rejected(self):
        self.output.symlink_to(self.original)
        with self.assertRaisesRegex(ValueError, 'aliases outside|linked destination'): self.imported()
        self.output.unlink()
        alias = self.cache / 'escaped'; alias.symlink_to(self.folder, target_is_directory=True)
        self.kwargs['output_pickle'] = alias / 'other.pkl'
        with self.assertRaisesRegex(ValueError, 'aliases outside|linked destination'): self.imported()

    def test_outside_root_destination_is_rejected(self):
        self.kwargs['output_pickle'] = self.folder / 'outside.pkl'
        with self.assertRaisesRegex(ValueError, 'leaves configured'): self.imported()

    def test_linked_receipt_preserves_its_target(self):
        self.ready.symlink_to(self.run_record)
        before = self.run_record.read_bytes()
        with self.assertRaisesRegex(ValueError, 'aliases outside|linked destination'): self.imported()
        self.assertEqual(self.run_record.read_bytes(), before)

    def test_output_change_before_receipt_does_not_close_publication(self):
        original_json = m.write_json
        def writer(path, value, **kwargs):
            if Path(path) == self.ready:
                self.output.write_bytes(b'changed after payload publication')
            return original_json(path, value, **kwargs)
        with patch.object(m, 'write_json', writer):
            with self.assertRaisesRegex(ValueError, 'changed'): self.imported()
        self.assertTrue(self.pending.exists())
        with self.assertRaisesRegex(ValueError, 'pending'): self.check()

    def test_changed_completed_stage_is_preserved(self):
        with self.interrupted('payload'):
            with self.assertRaises(OSError): self.imported()
        stage = Path(m.read_json(self.pending)['stage']); stage.write_bytes(b'corrupt stage')
        with self.assertRaisesRegex(ValueError, 'corrupted'): self.imported()
        self.assertEqual(stage.read_bytes(), b'corrupt stage')

    def test_interrupted_after_receipt_recovers_before_reuse(self):
        original_unlink = Path.unlink
        def unlink(path, *args, **kwargs):
            if path == self.pending:
                raise OSError('SYNTHETIC interruption before journal close')
            return original_unlink(path, *args, **kwargs)
        with patch.object(Path, 'unlink', unlink):
            with self.assertRaisesRegex(OSError, 'journal close'): self.imported()
        self.assertTrue(self.ready.exists())
        with self.assertRaisesRegex(ValueError, 'pending'): self.check()
        before = m.entry(self.output, hashed=True)
        self.imported()
        self.assertEqual(m.entry(self.output, hashed=True), before)
        self.assertTrue(self.check()['reused'])

    def test_installer_reuses_original_without_builder_or_yaml_write(self):
        self.imported(); result, events = self.shell()
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertEqual(events, '')
        self.assertIn('completed original', result.stderr)

    def test_installer_missing_original_receipt_uses_existing_builder_route(self):
        result, events = self.shell()
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertIn('BUILDER:', events); self.assertIn('UPDATE:', events)

    def test_installer_invalid_original_receipt_fails_without_builder(self):
        self.imported(); self.output.write_bytes(b'changed')
        result, events = self.shell()
        self.assertNotEqual(result.returncode, 0)
        self.assertEqual(events, '')

    def test_installer_pending_original_import_fails_without_builder(self):
        with self.interrupted('receipt'):
            with self.assertRaises(OSError): self.imported()
        result, events = self.shell()
        self.assertNotEqual(result.returncode, 0)
        self.assertEqual(events, '')

    def test_installer_preserves_explicit_external_legacy_reuse(self):
        outside = self.folder / 'external.pkl'; outside.write_bytes(b'legacy explicit selection')
        future = max(self.vcf.stat().st_mtime_ns, (ROOT / 'scripts/stat_protein_domain_amscores.py').stat().st_mtime_ns) + 1000000000
        os.utime(outside, ns=(future, future))
        result, events = self.shell(output=outside)
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertEqual(events, '')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--fixture-root', type=Path, required=True)
    parser.add_argument('--receipt', type=Path)
    args = parser.parse_args()
    global FIXTURE_ROOT
    FIXTURE_ROOT = args.fixture_root.absolute()
    FIXTURE_ROOT.mkdir(parents=True, exist_ok=False)
    result = unittest.TextTestRunner(verbosity=2).run(unittest.defaultTestLoader.loadTestsFromTestCase(OriginalCache))
    value = dict(schema='priva-original-domain-publication-fixtures-v1', checkout=str(ROOT),
                 tests_run=result.testsRun, failures=len(result.failures), errors=len(result.errors),
                 successful=result.wasSuccessful(), helper_sha256=m.digest(HELPER),
                 tests_sha256=m.digest(__file__), installer_sha256=m.digest(ROOT / 'scripts/install_utils.sh'),
                 fixture_root=str(FIXTURE_ROOT), finished_utc=datetime.now(timezone.utc).isoformat(),
                 scope='synthetic publication/routing only; no real cache load, producer run, or installation acceptance')
    if args.receipt:
        with args.receipt.open('x') as stream:
            json.dump(value, stream, indent=2); stream.write('\n')
    print(json.dumps(value, sort_keys=True))
    return 0 if result.wasSuccessful() else 1


if __name__ == '__main__':
    raise SystemExit(main())
