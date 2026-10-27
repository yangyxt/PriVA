"""Portable source/codec/process controls; real native parity is a separate run.

No model, cache or family input is used. Temporary outputs follow TMPDIR.
"""
import ast
import hashlib
import importlib
import io
import os
from pathlib import Path
import stat
import sys
import tempfile
from types import SimpleNamespace
import unittest
from unittest.mock import patch

import pysam
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'scripts'))
import pedigree_filtration_per_fam as reference
import priva_rust_pedigree_filtration as bridge


def fake_executable(root, mode='ok'):
    """Real child/pipe fixture, no scientific decision implementation."""
    path = root / ('fake-' + mode)
    source = '''#!PYTHON
import csv,io,os,struct,sys,time
MODE=MODE_LITERAL
def read():
 h=sys.stdin.buffer.read(8)
 if not h:return None
 n,=struct.unpack('>Q',h)
 rows=list(csv.reader(io.StringIO(sys.stdin.buffer.read(n).decode())))
 return dict(rows[1:])
def write(fields):
 s=io.StringIO(newline='');w=csv.writer(s,lineterminator='\\n');w.writerow(['name','value']);w.writerows(fields)
 b=s.getvalue().encode();sys.stdout.buffer.write(struct.pack('>Q',len(b))+b);sys.stdout.buffer.flush()
while True:
 q=read()
 if q is None:break
 if MODE=='hang':time.sleep(60)
 if MODE=='short':sys.stdout.buffer.write(b'bad');sys.stdout.buffer.flush();break
 seq=q['sequence']
 if MODE=='ordinal':seq='1.0'
 fields=[('status','ok'),('sequence',seq),('value','priva-pedigree-v1' if q['op']=='hello' else 'closed')]
 if MODE=='duplicate':fields.append(('value','duplicate'))
 if MODE=='handshake':fields[-1]=('value','different-profile')
 if q['op']=='genotype':fields=[('status','source-error'),('sequence',seq),('class','TypeError'),('message','fixture source error')]
 write(fields)
 if q['op']=='quit':
  if MODE=='extra':sys.stdout.buffer.write(b'extra');sys.stdout.buffer.flush()
  if MODE=='nonzero':sys.exit(37)
  break
'''.replace('PYTHON', sys.executable, 1).replace('MODE_LITERAL', repr(mode), 1)
    path.write_text(source)
    path.chmod(0o700)
    return path


class PedigreePortableControls(unittest.TestCase):
    def tearDown(self):
        bridge.close_sessions()

    def test_twelve_original_bodies_remain_admitted(self):
        bridge._admit(vars(reference))
        source = (ROOT / 'scripts/pedigree_filtration_per_fam.py').read_text()
        observed = {}
        for node in ast.parse(source).body:
            if isinstance(node, ast.FunctionDef):
                observed[node.name] = hashlib.sha256(ast.get_source_segment(source, node).encode()).hexdigest()
            elif isinstance(node, ast.ClassDef) and node.name == 'VariantSource':
                for method in node.body:
                    if isinstance(method, ast.FunctionDef):
                        observed['VariantSource.' + method.name] = hashlib.sha256(ast.get_source_segment(source, method).encode()).hexdigest()
        self.assertEqual(observed, bridge._EXPECTED_BODIES)
        self.assertEqual(len(observed), 12)

    def test_rebound_callable_is_not_silently_accepted(self):
        with patch.object(reference, 'is_variant_carrier', lambda value: True):
            with self.assertRaises(bridge.Unsupported):
                bridge._admit(vars(reference))

    def test_class_hooks_and_logger_callbacks_are_not_silently_skipped(self):
        with patch.object(reference.VariantSource, '__new__', staticmethod(lambda cls, *args: object.__new__(cls))):
            with self.assertRaisesRegex(bridge.Unsupported, 'class hooks'):
                bridge._admit(vars(reference))
        with patch.object(reference.logger, 'info', lambda *args: None):
            with self.assertRaisesRegex(bridge.Unsupported, 'logger instance'):
                bridge._admit(vars(reference))

    def test_direct_unsupported_source_does_not_start_a_child(self):
        with patch.object(reference, 'is_variant_carrier', lambda value: 'original-callback'):
            with patch.object(bridge, 'NativeSession', side_effect=AssertionError('must not start')):
                self.assertEqual(bridge.is_variant_carrier((1,), binary='not-opened'), 'original-callback')

    def test_unsupported_whole_source_falls_back_before_native_or_IO(self):
        with patch.object(reference, 'pedigree_filter', return_value='original-result') as original:
            with patch.object(bridge, 'NativeSession', side_effect=AssertionError('must not start')):
                self.assertEqual(bridge.pedigree_filter('v', 'p', 'F', 'o', binary='not-opened'), 'original-result')
        original.assert_called_once_with('v', 'p', 'F', 'o')

    def test_list_array_and_string_subclass_family_targets_use_original_before_effects(self):
        class Text(str):
            pass
        for target in (['F'], np.asarray(['F']), Text('F'), np.str_('F')):
            with self.subTest(kind=type(target).__name__), patch.object(reference, 'pedigree_filter', return_value='old') as original:
                with patch.object(bridge, 'NativeSession', side_effect=AssertionError('must not start')):
                    self.assertEqual(bridge.pedigree_filter('v', 'p', target, 'o', binary='not-opened'), 'old')
                self.assertIs(original.call_args.args[2], target)

    def test_raw_sex_equality_uses_provided_GTs_without_rereading(self):
        class Text(str):
            pass
        class Client:
            def __init__(self):
                self.admitted_namespaces = {}
                self.operations = []
            def request(self, operation, **fields):
                self.operations.append(operation)
                if operation == 'chromosome':
                    return {'value': 'X'}
                if operation == 'genotype':
                    return {'value': 'true'}
                raise AssertionError('raw sex must not enter native routing/source')
        record = SimpleNamespace(chrom='X', samples={}, alts=('C',), pos=1, id='unit', ref='A')
        for sex in (np.str_('1'), Text('1')):
            with self.subTest(kind=type(sex).__name__):
                client = Client()
                instance = bridge.VariantSource(record, 'p', 'f', 'm', sex, session=client)
                with patch.object(instance, '_get_genotype', side_effect=AssertionError('provided GTs must not be reread')):
                    self.assertEqual(instance._determine_chrX_source((1,), None, (1,)), 'MATERNAL')
                    self.assertEqual(instance._determine_chrY_source((1,), (1,)), 'PATERNAL')
                self.assertNotIn('route', client.operations)
                self.assertNotIn('source', client.operations)

    def test_whole_source_getters_keep_None_parent_identity_and_reject_foreign_sample(self):
        class Client:
            def __init__(self, foreign=False):
                self.admitted_namespaces = {}
                self.next_sample = 0
                self.foreign = foreign
            def request(self, op, **fields):
                if op == 'filter-begin':
                    return {'action': 'constructor-chrom'}
                if op == 'chromosome':
                    return {'value': 'Y'}
                if self.next_sample < 3:
                    sample = ('p', '', '')[self.next_sample]
                    self.next_sample += 1
                    return {'action': 'source-gt', 'sample': 'foreign' if self.foreign else sample}
                return {'action': 'write', 'value': 'UNKNOWN'}
        context = bridge._Context('p', None, None, '1', {'p'})
        record = SimpleNamespace(chrom='Y', samples={'p': {'GT': (1,)}}, alts=('C',))
        old = bridge.VariantSource._get_genotype
        for foreign in (False, True):
            observed = []
            def getter(instance, sample):
                observed.append(sample)
                return old(instance, sample)
            with self.subTest(foreign=foreign), patch.object(bridge.VariantSource, '_get_genotype', getter):
                if foreign:
                    with self.assertRaisesRegex(bridge.NativeProtocolError, 'getter identity'):
                        bridge._filter_record(Client(True), record, 1, context, vars(reference))
                    self.assertEqual(observed, [])
                else:
                    self.assertEqual(bridge._filter_record(Client(), record, 1, context, vars(reference)), 'UNKNOWN')
                    self.assertEqual(observed, ['p', None, None])

    def test_actual_pandas_vector_and_string_scalar_targets_retain_the_matching_family(self):
        class Text(str):
            pass
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            ped, vcf = root / 'small.ped', root / 'small.vcf'
            ped.write_text('F\tp\t0\t0\t1\t2\n')
            vcf.write_text('##fileformat=VCFv4.2\n##contig=<ID=1,length=10>\n'
                '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n'
                '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tp\n'
                '1\t1\tr1\tA\tC\t.\tPASS\t.\tGT\t0/1\n')
            for index, target in enumerate((['F'], np.asarray(['F']), Text('F'), np.str_('F'))):
                stats = {}
                destination = root / f'output-{index}.vcf'
                with self.subTest(kind=type(target).__name__):
                    with patch.object(bridge, 'NativeSession', side_effect=AssertionError('unsupported target must not start native')):
                        bridge.pedigree_filter(str(vcf), str(ped), target, str(destination), binary='not-opened', statistics=stats)
                    with pysam.VariantFile(destination) as result:
                        self.assertEqual([record.id for record in result], ['r1'])
                    self.assertEqual((stats['backend'], stats['native_processes']), ('python', 0))
                    self.assertTrue(stats['fallback_reason'])

    def test_default_whole_route_does_not_open_a_binary(self):
        with patch.object(reference, 'pedigree_filter', return_value=None) as original:
            with patch.object(bridge, 'NativeSession', side_effect=AssertionError('must not start')):
                bridge.pedigree_filter('v', 'p', 'F', 'o')
        original.assert_called_once()

    def test_typed_codec_keeps_none_empty_bool_int_and_string_distinct(self):
        values = [None, (), [], (None,), (True,), (1,), ('1',), 1, '0/1']
        encoded = [bridge._encode_gt(value) for value in values]
        self.assertEqual(len({tuple(row.items()) for row in encoded}), len(values))
        with self.assertRaises(bridge.Unsupported):
            bridge._encode_gt(([1],))

    def test_csv_round_trip_preserves_unicode_tabs_and_newlines(self):
        fields = dict(a='metadata Ω', b='one\ttwo\nthree', c='', d='1.0')
        framed = bridge._frame(fields)
        self.assertEqual(bridge._parse_frame(framed[8:]), fields)

    def test_csv_duplicate_empty_and_invalid_UTF8_are_rejected(self):
        for payload in (b'name,value\na,1\na,2\n', b'name,value\n,1\n', b'wrong,header\na,1\n', b'name,value\na,\xff\n'):
            with self.subTest(payload=payload), self.assertRaises(bridge.NativeProtocolError):
                bridge._parse_frame(payload)

    def test_real_child_has_one_start_and_an_explicit_close(self):
        with tempfile.TemporaryDirectory() as directory:
            with bridge.NativeSession(fake_executable(Path(directory))) as session:
                self.assertEqual(session.starts, 1)
                self.assertEqual(session.sequence, 1)
            self.assertTrue(session.closed)
            self.assertEqual(session.process.returncode, 0)

    def test_source_exception_does_not_become_protocol_failure_or_retry(self):
        with tempfile.TemporaryDirectory() as directory:
            with bridge.NativeSession(fake_executable(Path(directory))) as session:
                with self.assertRaisesRegex(TypeError, 'fixture source error'):
                    session.request('genotype')
                self.assertEqual(session.sequence, 2)
                self.assertEqual(session.starts, 1)

    def test_bad_handshake_short_duplicate_and_noncanonical_ordinal_fail(self):
        for mode in ('handshake', 'short', 'duplicate', 'ordinal'):
            with self.subTest(mode=mode), tempfile.TemporaryDirectory() as directory:
                with self.assertRaises(bridge.NativeProtocolError):
                    bridge.NativeSession(fake_executable(Path(directory), mode), timeout=0.5)

    def test_hanging_child_is_bounded_and_not_restarted(self):
        with tempfile.TemporaryDirectory() as directory:
            with self.assertRaisesRegex(bridge.NativeProtocolError, 'timed out'):
                bridge.NativeSession(fake_executable(Path(directory), 'hang'), timeout=0.05)

    def test_nonzero_or_extra_close_output_is_a_failure(self):
        for mode in ('nonzero', 'extra'):
            with self.subTest(mode=mode), tempfile.TemporaryDirectory() as directory:
                session = bridge.NativeSession(fake_executable(Path(directory), mode))
                with self.assertRaises(bridge.NativeProtocolError):
                    session.close()
                self.assertTrue(session.closed)

    def test_changed_executable_is_not_warm_reuse(self):
        with tempfile.TemporaryDirectory() as directory:
            path = fake_executable(Path(directory))
            session = bridge.NativeSession(path)
            path.write_text(path.read_text() + '\n# changed\n')
            with self.assertRaisesRegex(bridge.NativeProtocolError, 'executable changed'):
                session.close()

    def test_short_request_write_is_not_success(self):
        with tempfile.TemporaryDirectory() as directory:
            path = fake_executable(Path(directory))
            with patch.object(bridge.os, 'write', return_value=0):
                with self.assertRaisesRegex(bridge.NativeProtocolError, 'short request write'):
                    bridge.NativeSession(path, timeout=0.1)

    def test_actual_pysam_unknown_truthy_sample_raises_but_falsy_identifier_skips(self):
        header = pysam.VariantHeader()
        header.contigs.add('1')
        header.formats.add('GT', 1, 'String', 'Genotype')
        header.add_sample('p')
        record = header.new_record(contig='1', start=0, alleles=('A', 'C'))
        instance = bridge.VariantSource(record, 'p', None, None, '1')
        for sample in (None, ''):
            self.assertIsNone(instance._get_genotype(sample))
        with self.assertRaisesRegex(KeyError, 'invalid sample name: absent'):
            instance._get_genotype('absent')

    def test_getter_keeps_proxy_false_membership_separate(self):
        record = SimpleNamespace(chrom='1', samples={'p': {'GT': (1, None)}})
        instance = bridge.VariantSource(record, 'p', None, None, '1')
        self.assertIsNone(instance._get_genotype('absent'))
        self.assertEqual(instance._get_genotype('p'), (1, None))

    def test_invalid_binary_preserves_plain_file_and_rejects_empty_or_nonexecutable(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / 'ordinary'
            path.write_bytes(b'foreign')
            before = path.stat()
            with self.assertRaises(bridge.NativeProtocolError):
                bridge.NativeSession(path)
            self.assertEqual(path.read_bytes(), b'foreign')
            self.assertEqual(path.stat().st_ino, before.st_ino)
            path.write_bytes(b'')
            path.chmod(0o700)
            with self.assertRaises(bridge.NativeProtocolError):
                bridge.NativeSession(path)


if __name__ == '__main__':
    unittest.main()
