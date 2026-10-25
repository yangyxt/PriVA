"""Focused acquisition boundaries; all HTTP and payloads are synthetic."""
from contextlib import redirect_stderr, redirect_stdout
import hashlib
import importlib.util
import io
import json
import os
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch
from urllib.error import URLError

SOURCE = Path(__file__).resolve().parents[1] / 'scripts/cadd_refseq_resources.py'
spec = importlib.util.spec_from_file_location('cadd_refseq_under_test', SOURCE)
M = importlib.util.module_from_spec(spec)
spec.loader.exec_module(M)


class Response:
    def __init__(self, url, data=b'', status=200, headers=None, interrupt=False):
        self.url, self.data, self.status = url, data, status
        self.headers = {'Content-Length': str(len(data)), 'Last-Modified': 'Wed, 06 Aug 2025 00:00:00 GMT'}
        self.headers.update(headers or {})
        self.position = 0
        self.interrupt = interrupt

    def geturl(self): return self.url
    def __enter__(self): return self
    def __exit__(self, *args): return False

    def read(self, size=-1):
        if self.interrupt and self.position:
            raise URLError('synthetic interrupted transfer')
        if self.interrupt: size = 3
        stop = len(self.data) if size < 0 else self.position + size
        result = self.data[self.position:stop]; self.position += len(result)
        return result


class Publisher:
    def __init__(self, files, payloads):
        self.files, self.payloads = files, payloads
        self.calls = []
        self.interrupt = set()
        self.ignore_range = False
        self.wrong_range = False
        self.wrong_payload = False
        self.payload_encoding = None
        self.payload_last_modified = None
        self.manifest_encoding = None
        self.manifest_extra = b''
        self.manifest_wrong = False
        self.head_wrong = False
        self.redirect_head = False

    def __call__(self, request, timeout=60):
        self.calls.append(request)
        url = request.full_url
        if url == M.BASE + 'md5checksums.txt':
            data = ''.join(f'{md5}  ./{filename}\n' for filename, _, md5 in self.files.values()).encode()
            if self.manifest_wrong: data = data.replace(self.files['gff'][2].encode(), b'0' * 32)
            data += self.manifest_extra
            return Response(url, data, headers={'Content-Encoding': self.manifest_encoding})
        name = next(k for k, value in self.files.items() if url == M.BASE + value[0])
        data = self.payloads[name]
        if request.get_method() == 'HEAD':
            return Response(url + ('?other' if self.redirect_head else ''), headers={
                'Content-Length': str(len(data) + int(self.head_wrong))})
        byte_range = request.get_header('Range')
        headers = {'Content-Encoding': self.payload_encoding}
        if self.payload_last_modified: headers['Last-Modified'] = self.payload_last_modified
        status = 200
        if byte_range and not self.ignore_range:
            start = int(byte_range.removeprefix('bytes=').removesuffix('-'))
            headers['Content-Range'] = f'bytes {start}-{len(data) - 1}/{len(data)}'
            if self.wrong_range: headers['Content-Range'] = 'bytes 0-1/2'
            data = data[start:]; status = 206
        if self.wrong_payload: data = b'X' * len(data)
        interrupt = name in self.interrupt
        self.interrupt.discard(name)
        return Response(url, data, status, headers, interrupt)


class RefSeqAcquisitionTests(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory(prefix='priva-refseq-test-', dir=os.environ.get('TMPDIR'))
        self.addCleanup(self.temporary.cleanup)
        self.root = Path(self.temporary.name)
        self.cache = self.root / 'cache'
        self.info = self.cache / 'VEP_cache/homo_sapiens_merged/116_GRCh38/info.txt'
        self.info.parent.mkdir(parents=True)
        self.info.write_text('\n'.join([
            'assembly\tGRCh38', 'source_assembly\tGRCh38.p14', 'species\thomo_sapiens',
            'source_refseq\t' + M.RELEASE + ' - ' + M.PREFIX + '_genomic.gff',
            'bam\t' + M.PREFIX + '_knownrefseq_alns.bam', 'source_gencode\tGENCODE 50', '']))
        self.config = dict(cache_root=str(self.cache), assembly='hg38', vep_installed_version='116')
        self.payloads = {name: (name + ':0123456789').encode() for name in M.FILES}
        self.files = {name: (value[0], len(self.payloads[name]), hashlib.md5(self.payloads[name]).hexdigest())
                      for name, value in M.FILES.items()}
        self.pins_patch = patch.object(M, 'FILES', self.files)
        self.pins_patch.start(); self.addCleanup(self.pins_patch.stop)
        self.publisher = Publisher(self.files, self.payloads)
        self.network = patch.object(M, 'urlopen', self.publisher)
        self.network.start(); self.addCleanup(self.network.stop)

    def selected(self): return M.plan(self.config)
    def directory(self): return Path(self.selected()['directory'])
    def partial(self, name='gff'): return Path(self.selected()['resources'][name]['path'] + '.partial')

    def snapshot(self):
        return {str(p.relative_to(self.root)): (p.read_bytes(), M.identity(p))
                for p in self.root.rglob('*') if p.is_file() and not p.is_symlink()}

    def interrupt(self):
        self.publisher.interrupt.add('gff')
        with self.assertRaises(URLError): M.acquire(self.config)
        self.assertEqual(self.partial().read_bytes(), self.payloads['gff'][:3])

    def test_plan_is_local_readonly_and_does_not_create_candidate(self):
        before = self.snapshot(); result = self.selected()
        self.assertFalse(Path(result['directory']).exists())
        self.assertEqual(self.publisher.calls, [])
        self.assertEqual(before, self.snapshot())
        self.assertTrue(all(result[k] is False for k in M.FLAGS))

    def test_plan_requires_assembly_release_and_exact_cache_source(self):
        for updates in (dict(assembly='hg19'), dict(vep_installed_version='115'), dict(cache_root='relative')):
            with self.subTest(updates=updates), self.assertRaises(ValueError): M.plan({**self.config, **updates})
        self.info.write_text(self.info.read_text().replace('RS_2025_08', 'RS_2024_08'))
        with self.assertRaisesRegex(ValueError, 'no fallback'): self.selected()

    def test_duplicate_cache_metadata_fails(self):
        with self.info.open('a') as f: f.write('assembly\tGRCh38\n')
        with self.assertRaisesRegex(ValueError, 'Ambiguous'): self.selected()

    def test_linked_candidate_parent_is_refused(self):
        self.cache.joinpath('references').symlink_to(self.root, target_is_directory=True)
        with self.assertRaisesRegex(ValueError, 'linked candidate'): self.selected()

    def test_publisher_manifest_md5_length_and_encoding_guards(self):
        for attribute in ('manifest_wrong', 'head_wrong', 'redirect_head'):
            with self.subTest(attribute=attribute):
                setattr(self.publisher, attribute, True)
                with self.assertRaises(ValueError): M.publisher_pins(self.selected())
                setattr(self.publisher, attribute, False)
        self.publisher.manifest_encoding = 'gzip'
        with self.assertRaisesRegex(ValueError, 'Accept-Encoding'): M.publisher_pins(self.selected())

    def test_duplicate_publisher_row_is_refused(self):
        filename, _, md5 = self.files['gff']
        self.publisher.manifest_extra = f'{md5}  ./{filename}\n'.encode()
        with self.assertRaisesRegex(ValueError, 'ambiguous'): M.publisher_pins(self.selected())

    def test_acquire_records_raw_encoding_and_four_files_without_acceptance(self):
        result = M.acquire(self.config)
        self.assertEqual(set(result['resources']), set(self.files))
        for name, record in result['resources'].items():
            self.assertEqual(Path(record['path']).read_bytes(), self.payloads[name])
            self.assertEqual(record['publisher_md5'], self.files[name][2])
            self.assertEqual(record['sha256'], hashlib.sha256(self.payloads[name]).hexdigest())
            self.assertTrue(all(record[k] is False for k in M.FLAGS))
        source = M.read_json(self.directory() / 'source.ready.json')
        manifest = source['publisher_md5_manifest']
        self.assertEqual(manifest['http']['requested_accept_encoding'], 'identity')
        self.assertIsNone(manifest['http']['content_encoding'])
        self.assertEqual(manifest['sha256'], hashlib.sha256(Path(manifest['path']).read_bytes()).hexdigest())

    def test_check_is_local_readonly_and_hashes_completed_files(self):
        expected = M.acquire(self.config); before = self.snapshot()
        with patch.object(M, 'urlopen', side_effect=AssertionError('No HTTP in check')):
            self.assertEqual(M.check(self.config), expected)
        self.assertEqual(before, self.snapshot())

    def test_successful_reuse_does_not_make_network_requests(self):
        expected = M.acquire(self.config); before = self.snapshot(); self.publisher.calls.clear()
        self.assertEqual(M.acquire(self.config), expected)
        self.assertEqual(self.publisher.calls, [])
        # The existing lock file may change neither payloads nor receipts.
        self.assertEqual(before, self.snapshot())

    def test_interrupted_owned_partial_resumes_exact_range(self):
        self.interrupt(); result = M.acquire(self.config)
        calls = [r for r in self.publisher.calls if r.get_header('Range')]
        self.assertEqual(len(calls), 1)
        self.assertEqual(calls[0].get_header('Range'), 'bytes=3-')
        self.assertEqual(calls[0].get_header('Accept-encoding'), 'identity')
        self.assertEqual(calls[0].get_header('If-unmodified-since'), 'Wed, 06 Aug 2025 00:00:00 GMT')
        self.assertEqual(Path(result['resources']['gff']['path']).read_bytes(), self.payloads['gff'])

    def test_ignored_or_wrong_resume_preserves_existing_partial(self):
        self.interrupt(); before = self.partial().read_bytes()
        for attribute in ('ignore_range', 'wrong_range'):
            with self.subTest(attribute=attribute):
                setattr(self.publisher, attribute, True)
                with self.assertRaisesRegex(ValueError, 'exact resume'): M.acquire(self.config)
                self.assertEqual(self.partial().read_bytes(), before)
                setattr(self.publisher, attribute, False)

    def test_partial_replacement_is_refused_before_append(self):
        self.interrupt(); partial = self.partial(); replacement = partial.with_suffix('.replacement')
        replacement.write_bytes(partial.read_bytes()); os.replace(replacement, partial)
        self.publisher.calls.clear()
        with self.assertRaisesRegex(ValueError, 'identity changed'): M.acquire(self.config)
        self.assertEqual(self.publisher.calls, [])

    def test_same_inode_prefix_edit_is_refused(self):
        self.interrupt(); partial = self.partial()
        with partial.open('r+b') as stream: stream.write(b'X')
        with self.assertRaises(ValueError): M.acquire(self.config)
        self.assertEqual(partial.read_bytes(), b'X' + self.payloads['gff'][1:3])

    def test_same_inode_crash_growth_checks_prior_prefix_and_resumes(self):
        self.interrupt(); partial = self.partial()
        with partial.open('ab') as stream: stream.write(self.payloads['gff'][3:5])
        result = M.acquire(self.config)
        self.assertEqual(Path(result['resources']['gff']['path']).read_bytes(), self.payloads['gff'])
        self.assertEqual([r.get_header('Range') for r in self.publisher.calls if r.get_header('Range')], ['bytes=5-'])

    def test_crash_growth_with_changed_prior_prefix_is_refused(self):
        self.interrupt(); partial = self.partial()
        with partial.open('r+b') as stream:
            stream.write(b'XXX'); stream.seek(0, 2); stream.write(b'YY')
        with self.assertRaisesRegex(ValueError, 'prefix changed'): M.acquire(self.config)

    def test_unowned_partial_and_final_are_never_overwritten(self):
        directory = self.directory(); directory.mkdir(parents=True)
        self.partial().write_bytes(b'keep partial')
        with self.assertRaisesRegex(ValueError, 'Unowned'): M.acquire(self.config)
        self.assertEqual(self.partial().read_bytes(), b'keep partial')
        self.partial().unlink()
        final = Path(self.selected()['resources']['gff']['path']); final.write_bytes(self.payloads['gff'])
        with self.assertRaisesRegex(ValueError, 'Unidentified'): M.acquire(self.config)
        self.assertEqual(final.read_bytes(), self.payloads['gff'])

    def test_source_change_during_partial_reuse_is_refused(self):
        self.interrupt(); self.info.write_text(self.info.read_text() + '# changed identity\n')
        with self.assertRaisesRegex(ValueError, 'source/cache identity'): M.acquire(self.config)
        self.assertEqual(self.partial().read_bytes(), self.payloads['gff'][:3])

    def test_payload_encoding_or_source_change_preserves_partial(self):
        self.interrupt(); before = self.partial().read_bytes()
        self.publisher.payload_encoding = 'gzip'
        with self.assertRaisesRegex(ValueError, 'encoding'): M.acquire(self.config)
        self.assertEqual(self.partial().read_bytes(), before)
        self.publisher.payload_encoding = None
        self.publisher.payload_last_modified = 'different'
        with self.assertRaisesRegex(ValueError, 'source identity'): M.acquire(self.config)
        self.assertEqual(self.partial().read_bytes(), before)

    def test_corrupt_payload_is_held_without_final_publication(self):
        self.publisher.wrong_payload = True
        with self.assertRaisesRegex(ValueError, 'checksum differs'): M.acquire(self.config)
        self.assertFalse(Path(self.selected()['resources']['gff']['path']).exists())
        self.assertEqual(self.partial().read_bytes(), b'X' * len(self.payloads['gff']))

    def test_publication_without_ready_receipt_recovers_only_owned_inode(self):
        original = M.publish_json
        def interrupted(path, value):
            if path.name == 'gff.ready.json': raise OSError('synthetic publication interruption')
            return original(path, value)
        with patch.object(M, 'publish_json', interrupted), self.assertRaises(OSError): M.acquire(self.config)
        path = Path(self.selected()['resources']['gff']['path']); before = M.identity(path)
        result = M.acquire(self.config)
        self.assertEqual(result['resources']['gff']['identity'], before)

    def test_changed_final_checksum_or_receipt_fails_readonly_check(self):
        result = M.acquire(self.config); path = Path(result['resources']['gff']['path'])
        with path.open('r+b') as stream: stream.write(b'X')
        before = self.snapshot()
        with self.assertRaisesRegex(ValueError, 'checksum differs'): M.check(self.config)
        self.assertEqual(before, self.snapshot())

    def test_managed_payload_symlink_is_refused(self):
        result = M.acquire(self.config); path = Path(result['resources']['gff']['path'])
        other = self.root / 'original'; path.rename(other); path.symlink_to(other)
        with self.assertRaisesRegex(ValueError, 'regular candidate'): M.check(self.config)
        self.assertEqual(other.read_bytes(), self.payloads['gff'])

    def test_cli_named_plan_and_missing_check_do_not_create_candidate(self):
        stdout, stderr = io.StringIO(), io.StringIO()
        with redirect_stdout(stdout), redirect_stderr(stderr):
            self.assertEqual(M.main(['plan', '--config-json', json.dumps(self.config)]), 0)
            self.assertEqual(M.main(['check', '--config-json', json.dumps(self.config)]), 1)
        self.assertEqual(json.loads(stdout.getvalue())['refseq_release'], M.RELEASE)
        self.assertFalse(self.directory().exists())


if __name__ == '__main__': unittest.main()
