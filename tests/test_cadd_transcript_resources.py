"""Actual tiny HTTP/gzip/filesystem boundaries; no public assets or inference."""
import copy
import gzip
from http.server import BaseHTTPRequestHandler, HTTPServer
import importlib.util
import json
import os
from pathlib import Path
import tempfile
import threading
import unittest
from unittest.mock import patch
from urllib.parse import urlsplit


ROOT = Path(__file__).resolve().parents[1]
SPEC = importlib.util.spec_from_file_location('selected_cadd_transcripts', ROOT / 'scripts/cadd_transcript_resources.py')
C = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(C)


def bsd(data):
    value = 0
    for byte in data:
        value = ((value >> 1) | ((value & 1) << 15)) + byte
        value &= 65535
    return [value, (len(data) + 1023) // 1024]


class Candidates(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory(dir=os.environ.get('PRIVA_TEST_TMPDIR', os.environ.get('TMPDIR')))
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)
        reference = self.root / 'selected.fa'; reference.write_bytes(b'>chr1\nAC\n')
        Path(str(reference) + '.fai').write_text('chr1\t2\t6\t2\t3\n')
        self.config = {'assembly': 'hg38', 'vep_installed_version': '116',
                       'cache_root': str(self.root / 'cache'), 'ref_genome': str(reference)}
        self.assets, self.requests, self.mode = {}, [], {}
        owner = self
        class Handler(BaseHTTPRequestHandler):
            def log_message(self, *args): pass
            def do_HEAD(self): self.serve(False)
            def do_GET(self): self.serve(True)
            def serve(self, body):
                owner.requests.append((self.command, self.path, self.headers.get('Range')))
                if self.path == owner.mode.get('fail'):
                    self.send_error(503); return
                if self.path.endswith('/CHECKSUMS'):
                    prefix = self.path[:-len('CHECKSUMS')]
                    data = ''.join('%d %d %s\n' % (*bsd(value), Path(path).name)
                                   for path, value in owner.assets.items() if str(Path(path).parent) + '/' == prefix).encode()
                else:
                    data = owner.assets.get(self.path)
                    if data is None: self.send_error(404); return
                offset = 0
                if self.headers.get('Range') and not owner.mode.get('ignore_range'):
                    offset = int(self.headers['Range'].removeprefix('bytes=').removesuffix('-'))
                    self.send_response(206)
                    self.send_header('Content-Range', 'bytes %d-%d/%d' % (offset, len(data) - 1, len(data)))
                else: self.send_response(200)
                self.send_header('Content-Length', str(len(data) - offset)); self.send_header('ETag', '"tiny-fixed"')
                self.end_headers()
                if body:
                    try: self.wfile.write(data[offset:])
                    except (BrokenPipeError, ConnectionResetError): pass
        self.server = HTTPServer(('127.0.0.1', 0), Handler)
        self.thread = threading.Thread(target=self.server.serve_forever, daemon=True); self.thread.start()
        self.addCleanup(self.close_server)
        origin = 'http://127.0.0.1:%d/pub' % self.server.server_port
        p = patch.object(C, 'ENSEMBL', origin); p.start(); self.addCleanup(p.stop)
        self.seed()

    def close_server(self):
        self.server.shutdown(); self.thread.join(); self.server.server_close()

    def bodies(self, assembly='GRCh38'):
        return {'peptide': ('>ENSP1.1 pep chromosome:%s:1:1:3:1 gene:ENSG1.1 transcript:ENST1.1\nMA\n' % assembly).encode(),
                'gtf': ('#!genome-build %s.p14\n#!genome-version %s\n1\tensembl\texon\t1\t3\t.\t+\t.\tgene_id "G"; transcript_id "T";\n' % (assembly, assembly)).encode()}

    def seed(self):
        self.selected = C.plan(self.config)
        for name, resource in self.selected['resources'].items():
            self.assets[urlsplit(resource['url']).path] = gzip.compress(self.bodies(self.selected['assembly'])[name], mtime=0)

    def raw_path(self, component='peptide'): return Path(self.selected['resources'][component]['path'])
    def state(self, name): return Path(self.selected['directory']) / name
    def write_json(self, path, value):
        path.parent.mkdir(parents=True, exist_ok=True); path.write_text(json.dumps(value, sort_keys=True))

    def stage_partial(self, count=7):
        resource = self.selected['resources']['peptide']
        source = C.source_pin(resource)
        self.state('peptide.download.json').parent.mkdir(parents=True)
        C.publish_json(self.state('peptide.download.json'), source)
        data = self.assets[urlsplit(resource['url']).path]
        Path(str(self.raw_path()) + '.partial').write_bytes(data[:count])
        return data[:count]

    def legacy_receipts(self, component='peptide'):
        resource = self.selected['resources'][component]
        historical = self.root / 'historical'; historical.mkdir(exist_ok=True)
        archive = historical / Path(resource['path']).name
        archive.write_bytes(self.assets[urlsplit(resource['url']).path])
        view = historical / Path(resource['view']).name; view.write_bytes(self.bodies(self.selected['assembly'])[component])
        raw_path, view_path = historical / (component + '.raw.json'), historical / (component + '.view.json')
        raw = {'assembly': self.selected['assembly'], 'component': 'ensembl116_' + component,
               'stage': 'raw_acquisition', 'status': 'source_acquired_not_scientifically_qualified',
               'scientific_qualified': False, 'cadd_enabled': False,
               'files': [{'path': str(archive), 'identity': C.R.identity(archive), 'sha256': C.R.digest(archive)}],
               'source': {'url': resource['url'], 'publisher_checksums': resource['checksums_url'],
                          'publisher_checksum_manifest_sha256': 'a' * 64,
                          'publisher_bsd_sum': bsd(archive.read_bytes()), 'bytes': archive.stat().st_size}}
        self.write_json(raw_path, raw)
        prepared = {'assembly': self.selected['assembly'], 'component': 'ensembl116_' + component,
                    'stage': 'uncompressed_view', 'status': 'view_prepared_not_scientifically_qualified',
                    'scientific_qualified': False, 'cadd_enabled': False, 'path': str(view),
                    'identity': C.R.identity(view), 'sha256': C.R.digest(view),
                    'raw_receipt': str(raw_path), 'raw_receipt_sha256': C.R.digest(raw_path),
                    'operation': 'gzip_decompression_preserving_all_original_bytes_and_record_order'}
        self.write_json(view_path, prepared)
        return raw_path, C.R.digest(raw_path), view_path, C.R.digest(view_path)

    def test_grch38_real_http_install_produces_exact_plain_views_and_separate_candidate_manifest(self):
        original = copy.deepcopy(self.config)
        result = C.install(self.config)
        record = C.read_json(result['path'])
        self.assertEqual(self.config, original)
        self.assertEqual(set(record['resources']), {'peptide', 'gtf'})
        for name, resource in self.selected['resources'].items():
            self.assertEqual(Path(resource['view']).read_bytes(), self.bodies()[name])
        self.assertEqual(record['selected_reference']['path'], self.config['ref_genome'])
        self.assertEqual(record['selected_reference']['identity'], dict(zip(
            ('st_ino', 'st_size', 'st_mtime_ns', 'st_ctime_ns'), C.R.identity(Path(self.config['ref_genome']))[1:])))
        self.assertIn('node_local_st_dev_excluded', record['selected_reference']['identity_policy'])
        self.assertTrue(all(record[key] is False for key in C.FLAGS))
        self.assertEqual(result['environment_candidates']['CADD_MMSPLICE_GTF'], self.selected['resources']['gtf']['view'])

    def test_grch37_uses_requested_release_endpoint_and_retains_annotation_base75_without_gtf(self):
        self.config.update(assembly='hg19'); self.seed()
        result = C.install(self.config)
        self.assertEqual(result['vep_release'], '116'); self.assertEqual(result['annotation_base_release'], '75')
        self.assertNotIn('CADD_MMSPLICE_GTF', result['environment_candidates'])
        self.assertIn('/grch37/release-116/', self.selected['resources']['peptide']['url'])
        self.assertFalse(any('/gtf/' in path for _, path, _ in self.requests))

    def test_requested_release_is_not_a_fixed116_or_fallback(self):
        self.config['vep_installed_version'] = 113; self.seed()
        result = C.install(self.config)
        self.assertEqual(result['vep_release'], '113')
        self.assertTrue(all('/release-113/' in path for _, path, _ in self.requests))

    def test_bad_release_or_root_is_rejected_without_network_or_files(self):
        for key, value in (('vep_installed_version', '116.3'), ('vep_installed_version', True),
                           ('assembly', 'unknown'), ('cache_root', ''), ('cache_root', '/')):
            config = dict(self.config, **{key: value})
            with self.subTest(key=key, value=value), self.assertRaises(ValueError): C.plan(config)
        self.assertEqual(self.requests, [])

    def test_completed_reuse_and_verify_only_are_offline_and_keep_payload_inodes(self):
        first = C.install(self.config)
        identities = {name: C.R.identity(Path(resource['view'])) for name, resource in self.selected['resources'].items()}
        with patch.object(C, 'urlopen', side_effect=AssertionError('No network on completed reuse')):
            self.assertEqual(C.install(self.config), first)
            self.assertEqual(C.install(self.config, verify_only=True), first)
        self.assertEqual(identities, {name: C.R.identity(Path(resource['view'])) for name, resource in self.selected['resources'].items()})

    def test_reference_device_changes_preserve_exact_candidate_offline_without_reference_hashing(self):
        first = C.install(self.config)
        reference = Path(self.config['ref_genome']); references = (reference, Path(str(reference) + '.fai'))
        original_identity, original_digest = C.R.identity, C.R.digest
        def digest(path, *args, **kwargs):
            self.assertNotIn(Path(path), references, 'Reference/index must not be scanned for receipt reuse')
            return original_digest(path, *args, **kwargs)
        for target in references:
            def identity(path):
                value = original_identity(path)
                if Path(path) == target: value[0] += 10000
                return value
            with self.subTest(target=target.name), \
                 patch.object(C.R, 'identity', side_effect=identity), \
                 patch.object(C.R, 'digest', side_effect=digest), \
                 patch.object(C, 'urlopen', side_effect=AssertionError('No network on completed reuse')):
                self.assertEqual(C.install(self.config, verify_only=True), first)

    def test_reference_other_identity_changes_refuse_prior_candidate_without_reference_hashing(self):
        first = C.install(self.config)
        reference = Path(self.config['ref_genome']); references = (reference, Path(str(reference) + '.fai'))
        original_identity, original_digest = C.R.identity, C.R.digest
        def digest(path, *args, **kwargs):
            self.assertNotIn(Path(path), references, 'Reference/index must not be scanned for receipt reuse')
            return original_digest(path, *args, **kwargs)
        for target in references:
            for field, label in enumerate(('st_ino', 'st_size', 'st_mtime_ns', 'st_ctime_ns'), 1):
                def identity(path):
                    value = original_identity(path)
                    if Path(path) == target: value[field] += 1
                    return value
                with self.subTest(target=target.name, field=label), \
                     patch.object(C.R, 'identity', side_effect=identity), \
                     patch.object(C.R, 'digest', side_effect=digest), \
                     patch.object(C, 'urlopen', side_effect=AssertionError('No network on verification')):
                    with self.assertRaises(FileNotFoundError): C.install(self.config, verify_only=True)
        self.assertEqual(C.install(self.config, verify_only=True), first)

    def test_verify_only_missing_resource_does_not_download(self):
        with patch.object(C, 'urlopen', side_effect=AssertionError('No acquisition')):
            with self.assertRaisesRegex(ValueError, 'not installed'): C.install(self.config, verify_only=True)

    def test_http_failure_is_not_complete_and_retry_recovers_owned_staging(self):
        self.mode['fail'] = urlsplit(self.selected['resources']['peptide']['url']).path
        with self.assertRaises(OSError): C.install(self.config)
        self.assertEqual(list(Path(self.selected['directory']).glob('candidate.*.ready.json')), [])
        self.mode.clear(); C.install(self.config)

    def test_exact_http_range_resumes_partial_without_replacing_prefix(self):
        self.stage_partial()
        C.install(self.config)
        self.assertTrue(any(method == 'GET' and offset == 'bytes=7-' for method, _, offset in self.requests))
        self.assertFalse(Path(str(self.raw_path()) + '.partial').exists())

    def test_ignored_range_keeps_partial_and_failure_propagates(self):
        prefix = self.stage_partial(); self.mode['ignore_range'] = True
        with self.assertRaisesRegex(ValueError, 'exact resume'): C.install(self.config)
        self.assertEqual(Path(str(self.raw_path()) + '.partial').read_bytes(), prefix)
        self.assertFalse(self.raw_path().exists())

    def test_changed_publisher_pin_rejects_staged_bytes(self):
        self.stage_partial(); self.config['vep_installed_version'] = '117'
        # Different release has a separate directory; it cannot adopt/resume116.
        with self.assertRaisesRegex(ValueError, 'Missing/ambiguous publisher checksum'): C.install(self.config)
        self.assertTrue(Path(str(self.raw_path()) + '.partial').exists())

    def test_unknown_existing_final_is_retained_without_network(self):
        self.raw_path().parent.mkdir(parents=True); self.raw_path().write_bytes(b'manual bytes')
        with patch.object(C, 'urlopen', side_effect=AssertionError('Do not acquire over originals')):
            with self.assertRaisesRegex(ValueError, 'Unidentified'): C.install(self.config)
        self.assertEqual(self.raw_path().read_bytes(), b'manual bytes')

    def test_changed_completed_view_is_not_replaced(self):
        C.install(self.config)
        view = Path(self.selected['resources']['peptide']['view']); view.write_bytes(b'changed')
        with self.assertRaisesRegex(ValueError, 'resource changed'): C.install(self.config)
        self.assertEqual(view.read_bytes(), b'changed')

    def test_wrong_assembly_and_truncated_gzip_never_publish_view_or_candidate_completion(self):
        resource = self.selected['resources']['peptide']; key = urlsplit(resource['url']).path
        self.assets[key] = gzip.compress(self.bodies('GRCh37')['peptide'], mtime=0)
        with self.assertRaisesRegex(ValueError, 'selected assembly'): C.install(self.config)
        self.assertFalse(Path(resource['view']).exists())
        self.assertEqual(list(Path(self.selected['directory']).glob('candidate.*.ready.json')), [])

    def test_truncated_gzip_propagates_and_retains_diagnostics(self):
        resource = self.selected['resources']['peptide']; key = urlsplit(resource['url']).path
        self.assets[key] = self.assets[key][:-4]
        with self.assertRaisesRegex(ValueError, 'gzip decompression failed'): C.install(self.config)
        self.assertFalse(Path(resource['view']).exists())
        self.assertTrue(list(Path(resource['view']).parent.glob('*.partial')))

    def test_linked_managed_parent_is_refused_without_touching_target(self):
        outside = self.root / 'outside'; outside.mkdir()
        cache = Path(self.config['cache_root']); cache.mkdir()
        (cache / 'references').symlink_to(outside)
        with self.assertRaisesRegex(ValueError, 'escapes|linked'): C.install(self.config)
        self.assertEqual(list(outside.iterdir()), [])

    def test_completed_payload_leaf_links_cannot_claim_managed_ownership(self):
        C.install(self.config)
        for value in ('path', 'view'):
            leaf = Path(self.selected['resources']['peptide'][value]); outside = self.root / ('external-' + leaf.name)
            leaf.rename(outside); leaf.symlink_to(outside)
            with self.subTest(value=value), self.assertRaisesRegex(ValueError, 'linked managed candidate leaf'):
                C.install(self.config)
            self.assertTrue(outside.is_file()); leaf.unlink(); outside.rename(leaf)

    def test_adoption_refuses_existing_external_leaf_without_touching_original(self):
        args = self.legacy_receipts()
        source = Path(C.read_json(args[0])['files'][0]['path']); identity = C.R.identity(source)
        self.raw_path().parent.mkdir(parents=True); self.raw_path().symlink_to(source)
        with self.assertRaisesRegex(ValueError, 'linked managed candidate leaf'):
            C.adopt(self.config, 'peptide', *args)
        self.assertEqual(C.R.identity(source), identity)

    def test_explicit_adoption_reuses_closed_raw_and_view_without_network_or_decompression(self):
        raw, raw_sha, view, view_sha = self.legacy_receipts()
        originals = C.read_json(raw)['files'][0]['path'], C.read_json(view)['path']
        identities = [C.R.identity(Path(path)) for path in originals]
        with patch.object(C, 'urlopen', side_effect=AssertionError('No download')), \
             patch.object(C.R, 'source_digest', side_effect=AssertionError('No duplicate derivation')):
            C.adopt(self.config, 'peptide', raw, raw_sha, view, view_sha)
            C.adopt(self.config, 'peptide', raw, raw_sha, view, view_sha)  # Retry-safe.
        for source, target in zip(originals, (self.raw_path(), Path(self.selected['resources']['peptide']['view']))):
            self.assertFalse(Path(source).samefile(target)); self.assertEqual(Path(source).read_bytes(), target.read_bytes())
        self.assertEqual(identities, [C.R.identity(Path(path)) for path in originals])
        self.assertEqual(C.R.digest(raw), raw_sha); self.assertEqual(C.R.digest(view), view_sha)

    def test_adopting_both_components_then_main_install_is_offline(self):
        for component in self.selected['resources']:
            C.adopt(self.config, component, *self.legacy_receipts(component))
        with patch.object(C, 'urlopen', side_effect=AssertionError('No duplicate downloads')):
            result = C.install(self.config)
        self.assertTrue(all(result[key] is False for key in C.FLAGS))

    def test_wrong_supplied_receipt_pin_or_release_is_rejected_before_linking(self):
        raw, raw_sha, view, view_sha = self.legacy_receipts()
        with self.assertRaisesRegex(ValueError, 'receipt checksum'): C.adopt(self.config, 'peptide', raw, '0' * 64, view, view_sha)
        self.config['vep_installed_version'] = '113'
        with self.assertRaisesRegex(ValueError, 'assembly/release/scope'): C.adopt(self.config, 'peptide', raw, raw_sha, view, view_sha)
        self.assertFalse(self.raw_path().exists())

    def test_modified_original_and_conflicting_managed_output_are_never_overwritten(self):
        raw, raw_sha, view, view_sha = self.legacy_receipts()
        self.raw_path().parent.mkdir(parents=True); self.raw_path().write_bytes(b'conflict')
        with self.assertRaisesRegex(ValueError, 'managed candidate differs'): C.adopt(self.config, 'peptide', raw, raw_sha, view, view_sha)
        self.assertEqual(self.raw_path().read_bytes(), b'conflict')
        Path(C.read_json(view)['path']).write_bytes(b'modified original')
        with self.assertRaisesRegex(ValueError, 'original changed'): C.adopt(self.config, 'peptide', raw, raw_sha, view, view_sha)

    def test_adoption_never_hardlinks_supplied_originals(self):
        args = self.legacy_receipts(); original_link = C.os.link
        def link(source, destination):
            if Path(source).parent.name == 'historical': raise AssertionError('Never hardlink supplied originals')
            return original_link(source, destination)
        with patch.object(C.os, 'link', side_effect=link): C.adopt(self.config, 'peptide', *args)
        self.assertEqual(self.raw_path().read_bytes(), Path(C.read_json(args[0])['files'][0]['path']).read_bytes())

    def test_source_content_change_during_copy_never_publishes_ready(self):
        args = self.legacy_receipts(); original_copy = C.shutil.copyfileobj
        def changed(input_stream, output, length):
            original_copy(input_stream, output, length)
            Path(input_stream.name).write_bytes(b'unapproved concurrent change')
        with patch.object(C.shutil, 'copyfileobj', side_effect=changed), self.assertRaisesRegex(ValueError, 'Original changed during copy'):
            C.adopt(self.config, 'peptide', *args)
        self.assertFalse(self.state('peptide.raw.ready.json').exists())


if __name__ == '__main__': unittest.main()
