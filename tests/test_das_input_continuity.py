"""Tiny provenance/installer controls, not scientific consumer-equivalence proof."""
import importlib.util
import json
import os
from pathlib import Path
import sys
import tempfile
from types import SimpleNamespace
import unittest
from unittest import mock

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT/'scripts'))
import das_input_continuity as C


def dump(path, value):
    path.write_text(json.dumps(value, sort_keys=True, indent=2) + '\n')


def certificate(folder, args, old, new):
    evidence = {}
    for role in ('plan', 'comparison', 'publication'):
        path = folder/(role+'.json')
        dump(path, {'fixture_only': True, 'role': role})
        evidence[role] = C.pin(path)
    sources = C.state.sources(args.source)
    value = dict(kind=C.PROOF_KIND, schema_version=1, status='completed', profile=args.profile,
        old_input=C.state.identities([old])[0], new_input=C.state.identities([new])[0],
        gates={name: True for name in ('full_record_census_equal', 'full_global_order_equal',
              'protected_annotation_fields_equal', 'headers_reviewed', 'publication_identity_exact')},
        consumer=dict(role='das_and_dependent_pm1', exact=True, differing_rows=0, ordered_rows=3,
            old_sha256='a'*64, new_sha256='a'*64,
            analysis_source_sha256=next(item['sha256'] for item in sources if item['role']=='analysis')),
        evidence=evidence)
    path = folder/'equivalence.json'
    dump(path, value)
    args.proof, args.proof_sha256 = str(path), C.pin(path)['sha256']
    return path, value


class InputContinuity(unittest.TestCase):
    def setUp(self):
        temporary = tempfile.TemporaryDirectory(prefix='das-input-continuity-',
            dir=os.environ.get('PRIVA_TEST_TMPDIR', os.environ.get('TMPDIR')))
        self.addCleanup(temporary.cleanup)
        self.root = Path(temporary.name)
        def file(name):
            path = self.root/name
            path.write_text('tiny fixture: '+name+'\n')
            return path
        self.old, self.new = file('old.vcf'), file('new.vcf')
        self.other, self.bundle, self.source = file('am.pkl'), file('bundle'), file('analysis.py')
        self.outputs = [file('out'+str(i)) for i in range(3)]
        self.receipt = self.root/'original.json'
        self.args = SimpleNamespace(profile='hg38', input=[str(self.other), str(self.old)],
            output=list(map(str, self.outputs)), source=['analysis='+str(self.source)],
            bundled=['curated='+str(self.bundle)], parameter=['threads=22', 'filter=true'],
            receipt=str(self.receipt), continuity=str(self.root/'continuity.json'),
            old_input=str(self.old), new_input=str(self.new))
        old = {'version': C.state.VERSION, **C.state.input_state(self.args),
               'outputs': C.state.identities(self.args.output)}
        dump(self.receipt, old)
        self.original = self.receipt.read_bytes()
        self.args.input[-1] = str(self.new)
        self.proof, self.value = certificate(self.root, self.args, self.old, self.new)

    def update_proof(self):
        dump(self.proof, self.value)
        self.args.proof_sha256 = C.pin(self.proof)['sha256']

    def reject_record(self):
        with self.assertRaises((OSError, ValueError, KeyError, TypeError)):
            C.record(self.args)
        self.assertFalse(Path(self.args.continuity).exists())
        self.assertEqual(self.receipt.read_bytes(), self.original)

    def test_record_check_and_idempotent_reuse_preserve_original_outputs(self):
        outputs = [path.read_bytes() for path in self.outputs]
        C.record(self.args)
        before = C.pin(self.args.continuity)
        C.check(self.args)
        C.record(self.args)
        self.assertEqual(C.pin(self.args.continuity), before)
        self.assertEqual(self.receipt.read_bytes(), self.original)
        self.assertEqual([path.read_bytes() for path in self.outputs], outputs)
        self.assertIs(json.loads(Path(self.args.continuity).read_text())['producer_execution'], False)

    def test_no_annotation_output_or_bundled_payload_read(self):
        protected = {self.old, self.new, self.other, self.bundle, *self.outputs}
        metadata = C.metadata
        def read(path):
            self.assertNotIn(Path(path), protected)
            return metadata(path)
        with mock.patch.object(C, 'metadata', side_effect=read):
            C.record(self.args)
            C.check(self.args)

    def test_each_required_gate_is_explicit_and_true(self):
        for name in self.value['gates']:
            with self.subTest(gate=name):
                self.value['gates'][name] = False
                self.update_proof()
                self.reject_record()
                self.value['gates'][name] = True

    def test_row_count_missingness_digest_and_source_fail_closed(self):
        for field, bad in [('exact', False), ('differing_rows', 1), ('differing_rows', False),
                           ('ordered_rows', 0), ('ordered_rows', True), ('new_sha256', 'b'*64),
                           ('analysis_source_sha256', 'c'*64), ('role', 'another-consumer')]:
            with self.subTest(field=field, value=bad):
                original = self.value['consumer'][field]
                self.value['consumer'][field] = bad
                self.update_proof()
                self.reject_record()
                self.value['consumer'][field] = original

    def test_proof_metadata_and_every_underlying_evidence_pin_are_bound(self):
        C.record(self.args)
        (self.root/'comparison.json').write_text('{"changed":true}\n')
        with self.assertRaisesRegex(ValueError, 'Pinned metadata changed'):
            C.check(self.args)

    def test_underlying_frozen_source_pin_does_not_require_json(self):
        source = self.root/'comparison-controller.py'
        source.write_text('# frozen fixture controller, not JSON\n')
        self.value['evidence']['controller'] = C.pin(source)
        self.update_proof()
        C.record(self.args)
        C.check(self.args)

    def test_changed_proof_without_new_explicit_checksum_is_rejected(self):
        self.proof.write_text(self.proof.read_text()+' ')
        self.reject_record()

    def test_incomplete_or_wrong_profile_proof_is_rejected(self):
        for name, bad in [('kind', 'other'), ('status', 'pending'),
                          ('profile', 'hg19'), ('schema_version', True)]:
            with self.subTest(field=name):
                old = self.value[name]
                self.value[name] = bad
                self.update_proof()
                self.reject_record()
                self.value[name] = old

    def test_original_receipt_and_nonclinvar_dependencies_stay_strict(self):
        for path in (self.other, self.source, self.outputs[0], self.old, self.new):
            with self.subTest(path=path.name):
                before = path.read_bytes(), path.stat()
                path.write_bytes(b'changed dependency\n')
                self.reject_record()
                path.write_bytes(before[0])
                os.utime(path, ns=(before[1].st_atime_ns, before[1].st_mtime_ns))

    def test_bundled_data_remain_existence_nonzero_only(self):
        C.record(self.args)
        self.bundle.write_bytes(b'different nonempty snapshot bytes')
        C.check(self.args)
        self.bundle.write_bytes(b'')
        with self.assertRaises(ValueError):
            C.check(self.args)

    def test_ordered_parameters_and_input_ordinal_cannot_change(self):
        self.args.parameter.reverse()
        self.reject_record()
        self.args.parameter.reverse()
        self.args.input.reverse()
        self.reject_record()

    def test_duplicate_new_input_is_rejected(self):
        self.args.input.append(str(self.new))
        self.reject_record()

    def test_receipt_aliases_and_links_are_rejected(self):
        self.args.continuity = str(self.outputs[0])
        with self.assertRaises(ValueError):
            C.record(self.args)
        self.assertTrue(self.outputs[0].exists())
        self.args.continuity = str(self.root/'continuity.json')
        self.proof.rename(self.root/'actual-proof.json')
        self.proof.symlink_to(self.root/'actual-proof.json')
        self.reject_record()

    def test_input_change_during_proof_validation_never_publishes(self):
        original = C.proof
        def changed(*args, **kwargs):
            result = original(*args, **kwargs)
            self.other.write_bytes(b'changed while checking evidence\n')
            return result
        with mock.patch.object(C, 'proof', side_effect=changed):
            self.reject_record()

    def test_evidence_change_after_its_first_scan_blocks_warm_reuse(self):
        C.record(self.args)
        original = C.proof
        def changed(*args, **kwargs):
            result = original(*args, **kwargs)
            (self.root/'comparison.json').write_text('{"changed_after_scan":true}\n')
            return result
        with mock.patch.object(C, 'proof', side_effect=changed), self.assertRaises(ValueError):
            C.check(self.args)

    def test_evidence_race_at_second_publication_check_never_publishes(self):
        original = C.proof
        calls = []
        def changed(*args, **kwargs):
            result = original(*args, **kwargs)
            calls.append(True)
            if len(calls) == 2:
                (self.root/'comparison.json').write_text('{"changed_at_publish":true}\n')
            return result
        with mock.patch.object(C, 'proof', side_effect=changed):
            self.reject_record()
        self.assertEqual(len(calls), 2)
        self.assertEqual(list(self.root.glob('continuity.json.*')), [])

    def test_publication_failure_preserves_old_receipt_and_removes_temporary(self):
        with mock.patch.object(C.os, 'link', side_effect=OSError('injected link failure')):
            self.reject_record()
        self.assertEqual(list(self.root.glob('continuity.json.*')), [])

    def test_unknown_receipt_fields_and_boolean_version_do_not_reuse(self):
        C.record(self.args)
        path = Path(self.args.continuity)
        value = json.loads(path.read_text())
        for updated in (dict(value, extra=True), dict(value, schema_version=True),
                        dict(value, producer_execution=0)):
            dump(path, updated)
            with self.assertRaises(ValueError):
                C.check(self.args)

    def test_real_shell_reuses_valid_continuity_without_scientific_writer(self):
        spec = importlib.util.spec_from_file_location('existing_das_shell_fixture',
            ROOT/'tests/test_das_install_completion.py')
        module = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(module)
        fixture = module.DASInstallCompletion()
        fixture.setUp()
        self.addCleanup(fixture.doCleanups)
        _, _, old = fixture.filter_inputs()
        fixture.assert_ok(fixture.run_install())
        original = fixture.receipt.read_bytes()
        completed = json.loads(original)
        new = fixture.folder/'new-clinvar.vcf.gz'
        new.write_bytes(b'new opaque ClinVar fixture\n')
        config = json.loads(fixture.config_path.read_text())
        config['clinvar_vcf'] = str(new)
        dump(fixture.config_path, config)
        roles = [('analysis','am_pick_intolerant_domains.py'), ('completion','das_install_state.py'),
                 ('mapping-io','interpro_mapping_io.py'), ('domain-mapping','protein_domain_mapping.py'),
                 ('statistics','stat_protein_domain_amscores.py'), ('domain-metadata','domain_metadata_cache.py')]
        bundle_paths = {'domain-normalization': ROOT/'data/InterPro/InterPro_domain_norm_map.json',
                        'functional-domains': ROOT/'data/InterPro/curated_InterPro_func_domains.tsv.gz',
                        'hcseeker-spots': fixture.folder/'bundled data/HCSeeker/HC_spots.hg38.tsv'}
        args = SimpleNamespace(profile=completed['profile'],
            input=[str(new) if row['path']==str(old) else row['path'] for row in completed['inputs']],
            output=[row['path'] for row in completed['outputs']],
            source=[role+'='+str(ROOT/'scripts'/name) for role,name in roles],
            bundled=[row['role']+'='+str(bundle_paths[row['role']]) for row in completed['bundled']],
            parameter=['='.join(pair) for pair in completed['parameters']],
            receipt=str(fixture.receipt), continuity=str(fixture.receipt)+'.clinvar-continuity.json',
            old_input=str(old), new_input=str(new))
        certificate(fixture.folder, args, old, new)
        C.record(args)
        fixture.assert_ok(fixture.run_install(mode='0'))
        self.assertEqual(len(fixture.calls()), 1)
        self.assertEqual(fixture.receipt.read_bytes(), original)
        self.assertEqual(json.loads(fixture.config_path.read_text())['clinvar_vcf'], str(new))


if __name__ == '__main__':
    unittest.main()
