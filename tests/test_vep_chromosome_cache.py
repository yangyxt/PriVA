"""Ordered-occurrence and live-prefix safety regressions; no VEP/model claim."""
import argparse
from contextlib import redirect_stdout
import gzip
import hashlib
import importlib.util
import io
import json
from pathlib import Path
import tempfile
import unittest
from unittest import mock


SPEC = importlib.util.spec_from_file_location(
    "vep_chromosome_cache", Path(__file__).parents[1] / "scripts/vep_chromosome_cache.py")
M = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(M)
HEADER = (b'##fileformat=VCFv4.2\n##contig=<ID=chr1,length=1000>\n'
          b'##contig=<ID=chr2,length=1000>\n##contig=<ID=NT_1.1,length=100>\n'
          b'##INFO=<ID=OLD,Number=1,Type=String,Description="source">\n'
          b'#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE\n')
ANNOTATED_HEADER = HEADER.replace(b'#CHROM', b'##INFO=<ID=CSQ,Number=.,Type=String,Description="Allele|Feature|HGNC_ID">\n#CHROM')
ROWS = [b'chr1\t10\tdup\tGG\tAA,C\t.\tPASS\tOLD=a\tGT\t./.\n',
        b'chr1\t10\tdup\tGG\tAA,C\t.\tPASS\tOLD=a\tGT\t./.\n',
        b'chr2\t20\tbnd\tA\tA]chr1:10]\t5\t.\t.\tGT\t0/1\n',
        b'NT_1.1\t2\textra\tC\tT\t.\tPASS\tOLD=other\tGT\t1\n']


def annotated(row, index):
    fields = row.rstrip(b'\n').split(b'\t')
    csq = b'CSQ=A|ENST' + str(index).encode() + b'|'
    fields[7] = csq if fields[7] == b'.' else fields[7] + b';' + csq
    return b'\t'.join(fields) + b'\n'


class PrefixTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)
        self.source = self.root / 'source.vcf.gz'
        with gzip.open(self.source, 'wb') as f:
            f.write(HEADER + b''.join(ROWS))
        self.serial = self.root / 'serial.vcf'
        self.serial.write_bytes(ANNOTATED_HEADER + b''.join(annotated(r, i) for i, r in enumerate(ROWS)))
        self.inventory = self.root / 'inventory.json'
        self.indexed = dict(input=M.identity(self.source), indexes=[],
                            header_sha256=hashlib.sha256(HEADER).hexdigest(),
                            contig_order=['chr1', 'chr2', 'NT_1.1'],
                            records_by_contig={'chr1': 2, 'chr2': 1, 'NT_1.1': 1}, total_records=4)
        self.write_inventory()

    def write_inventory(self):
        value = dict(kind='priva-vep-chromosome-inventory-v1',
                     input_inventory=self.indexed,
                     serial_output=dict(after=M.identity(self.serial),
                                        header_sha256=hashlib.sha256(ANNOTATED_HEADER).hexdigest()),
                     candidate_completed_contigs=['chr1', 'chr2'])
        self.inventory.write_text(json.dumps(value))

    def args(self):
        return argparse.Namespace(inventory=self.inventory,
                                  inventory_sha256=M.digest(self.inventory),
                                  output_dir=self.root/'preserved', cache_root=self.root,
                                  vep_source_root=getattr(self, 'vep_source_root', None))

    def preserve(self):
        with redirect_stdout(io.StringIO()):
            M.preserve_prefix(self.args())
        return json.loads((self.root/'preserved/preserved-prefix.json').read_text())

    def nonvariant_fixture(self, rows, emitted_indices, command=b'vep --offline --vcf --use_transcript_ref'):
        """Small independently pinned stand-ins exercise source/header gates.

        Production uses the four literal installed-module hashes in the helper;
        these fixtures do not claim to execute VEP or qualify an annotation.
        """
        self.vep_source_root = self.root/'vep-source'
        pins = {}
        for relative in M.VEP116_OMISSION_SOURCES:
            path = self.vep_source_root/relative
            path.parent.mkdir(parents=True, exist_ok=True)
            payload = ('pinned test source: ' + relative + '\n').encode()
            path.write_bytes(payload)
            pins[relative] = hashlib.sha256(payload).hexdigest()
        patcher = mock.patch.object(M, 'VEP116_OMISSION_SOURCES', pins)
        patcher.start()
        self.addCleanup(patcher.stop)
        header = ANNOTATED_HEADER.replace(
            b'#CHROM\t', b'##VEP="v116" API="v116"\n##VEP-command-line=\'' + command + b"'\n#CHROM\t")
        with gzip.open(self.source, 'wb') as f:
            f.write(HEADER + b''.join(rows))
        self.serial.write_bytes(header + b''.join(annotated(rows[i], i) for i in emitted_indices))
        counts = {}
        for row in rows:
            chrom = row.split(b'\t', 1)[0].decode()
            counts[chrom] = counts.get(chrom, 0) + 1
        self.indexed.update(input=M.identity(self.source), records_by_contig=counts,
                            total_records=len(rows))
        self.write_inventory()
        value = json.loads(self.inventory.read_text())
        value['serial_output']['header_sha256'] = hashlib.sha256(header).hexdigest()
        self.inventory.write_text(json.dumps(value))
        return header

    def test_exact_duplicates_multiallelic_breakends_missingness_and_samples(self):
        before = self.serial.read_bytes()
        proof = self.preserve()
        self.assertEqual(proof['validated_completed_contigs'], ['chr1', 'chr2'])
        self.assertEqual(proof['records_by_contig'], {'chr1': 2, 'chr2': 1})
        self.assertEqual((self.root/'preserved/contig-000.vep.vcf').read_bytes(),
                         ANNOTATED_HEADER + annotated(ROWS[0], 0) + annotated(ROWS[1], 1))
        self.assertEqual((self.root/'preserved/contig-001.vep.vcf').read_bytes(),
                         ANNOTATED_HEADER + annotated(ROWS[2], 2))
        self.assertEqual([(x['source_ordinal_start'],x['source_ordinal_end']) for x in proof['files']], [(1,2),(3,3)])
        self.assertEqual(self.serial.read_bytes(), before)
        self.assertFalse(proof['canonical_publication'])
        self.assertFalse(proof['scientific_qualification'])

    def test_allow_append_after_inventoried_extent(self):
        with self.serial.open('ab') as f:
            f.write(b'NT_1.1\t5\tuncertified-partial')
        proof = self.preserve()
        self.assertEqual(proof['records_by_contig']['chr1'], 2)
        self.assertTrue(proof['serial_writer_may_still_be_running'])

    def test_no_ready_proof_on_changed_allele(self):
        self.serial.write_bytes(self.serial.read_bytes().replace(b'GG\tAA,C', b'GG\tAT,C', 1))
        self.write_inventory()
        with self.assertRaisesRegex(ValueError, 'changed prefix occurrence'):
            self.preserve()
        self.assertFalse((self.root/'preserved/preserved-prefix.json').exists())

    def test_extra_non_vep_info_is_not_discarded(self):
        self.serial.write_bytes(self.serial.read_bytes().replace(b'OLD=a;CSQ=', b'OLD=a;NEW=1;CSQ=', 1))
        self.write_inventory()
        with self.assertRaisesRegex(ValueError, 'changed prefix occurrence'):
            self.preserve()

    def test_missing_info_marker_cannot_be_mixed_with_annotation(self):
        bad = annotated(ROWS[2],2).replace(b'\tCSQ=',b'\t.;CSQ=',1)
        with self.assertRaisesRegex(ValueError, 'INFO token mixed'):
            M.row_invariant(bad,True)

    def test_empty_info_tokens_are_not_normalized_away(self):
        for value in [b'OLD=a;;CSQ=A|ENST|',b';CSQ=A|ENST|',b'CSQ=A|ENST|;']:
            with self.subTest(info=value), self.assertRaisesRegex(ValueError, 'INFO token mixed'):
                fields=annotated(ROWS[2],2).split(b'\t'); fields[7]=value
                M.row_invariant(b'\t'.join(fields),True)

    def test_corrupted_destination_has_no_ready_proof(self):
        original_digest=M.digest
        def corrupt_copy(path):
            path=Path(path)
            if path.name.startswith('contig-'):
                with path.open('ab') as f: f.write(b'CORRUPT')
            return original_digest(path)
        with mock.patch.object(M,'digest',side_effect=corrupt_copy):
            with self.assertRaisesRegex(ValueError, 'copy checksum differs'):
                self.preserve()
        self.assertTrue((self.root/'preserved/contig-000.vep.vcf').exists())
        self.assertFalse((self.root/'preserved/preserved-prefix.json').exists())

    def test_partial_certified_record_fails_and_preserves_partial_files(self):
        self.serial.write_bytes(ANNOTATED_HEADER + annotated(ROWS[0],0) + b'chr1\t10')
        self.write_inventory()
        with self.assertRaisesRegex(ValueError, 'Incomplete VCF occurrence'):
            self.preserve()
        self.assertTrue((self.root/'preserved/contig-000.vep.vcf').exists())
        self.assertFalse((self.root/'preserved/preserved-prefix.json').exists())

    def test_missing_future_boundary_does_not_certify_prefix(self):
        self.serial.write_bytes(ANNOTATED_HEADER + b''.join(annotated(r,i) for i,r in enumerate(ROWS[:3])))
        self.write_inventory()
        with self.assertRaisesRegex(ValueError, 'Incomplete VCF occurrence'):
            self.preserve()

    def test_noncontiguous_source_contig_fails(self):
        with gzip.open(self.source,'wb') as f:
            f.write(HEADER + b''.join([ROWS[0],ROWS[2],ROWS[1],ROWS[3]]))
        self.indexed['input'] = M.identity(self.source)
        self.write_inventory()
        with self.assertRaises(ValueError):
            self.preserve()

    def test_changed_sample_header_fails(self):
        self.serial.write_bytes(self.serial.read_bytes().replace(b'\tSAMPLE\n',b'\tDIFFERENT\n',1))
        self.write_inventory()
        value = json.loads(self.inventory.read_text())
        value['serial_output']['header_sha256'] = hashlib.sha256(
            ANNOTATED_HEADER.replace(b'\tSAMPLE\n',b'\tDIFFERENT\n',1)).hexdigest()
        self.inventory.write_text(json.dumps(value))
        with self.assertRaisesRegex(ValueError, 'column/sample header'):
            self.preserve()

    def test_inventory_checksum_failure_has_no_output_directory(self):
        args = self.args(); args.inventory_sha256 = '0'*64
        with self.assertRaisesRegex(ValueError, 'Inventory checksum'):
            M.preserve_prefix(args)
        self.assertFalse(args.output_dir.exists())

    def test_source_identity_change_fails(self):
        with self.source.open('ab') as f:
            f.write(b'changed')
        with self.assertRaisesRegex(ValueError, 'Inventoried source changed'):
            self.preserve()

    def test_existing_destination_is_not_reused(self):
        (self.root/'preserved').mkdir()
        with self.assertRaisesRegex(ValueError, 'Fresh preservation directory'):
            self.preserve()

    def test_destination_must_stay_inside_cache_root(self):
        with self.assertRaisesRegex(ValueError, 'escapes cache_root'):
            M.owned(self.root.parent/'elsewhere',self.root)

    def test_inventory_is_explicitly_provisional_including_noncanonical_contigs(self):
        args = argparse.Namespace(input_vcf=self.source,serial_output=self.serial,
                                  output_json=self.root/'observed.json',cache_root=self.root,
                                  bcftools='unused',producer_job='1522950',cpus_per_chromosome=10,max_parallel=12)
        with mock.patch.object(M,'index_inventory',return_value=self.indexed), redirect_stdout(io.StringIO()):
            M.inventory(args)
        value=json.loads(args.output_json.read_text())
        self.assertEqual(value['candidate_completed_contigs'],['chr1','chr2'])
        self.assertEqual(value['active_contig'],'NT_1.1')
        self.assertEqual(value['validated_completed_contigs'],[])
        self.assertFalse(value['prefix_body_validated'])
        self.assertEqual(value['states']['chr1'],'observed_serial_prefix_NOT_VALIDATED')

    def test_pinned_nonvariant_omissions_keep_exact_duplicate_occurrences_and_source_ordinals(self):
        missing = b'chr1\t10\treference-only\tT\t.\t.\tPASS\tOLD=keep\tGT\t0/0\n'
        rows = [missing, ROWS[0], missing, missing, ROWS[1], missing, ROWS[2], ROWS[3]]
        header = self.nonvariant_fixture(rows, [1, 4, 6, 7])
        original = self.serial.read_bytes()
        proof = self.preserve()
        self.assertEqual(proof['records_by_contig'], {'chr1': 2, 'chr2': 1})
        self.assertEqual(proof['source_records_by_contig'], {'chr1': 6, 'chr2': 1})
        self.assertEqual(proof['omitted_records_by_contig'], {'chr1': 4, 'chr2': 0})
        ledger = proof['omitted_occurrences']
        path = Path(ledger['identity']['path'])
        omissions = [json.loads(line) for line in path.read_text().splitlines()]
        self.assertEqual([v['source_ordinal'] for v in omissions], [1, 3, 4, 6])
        self.assertEqual([v['source_contig_ordinal'] for v in omissions], [1, 3, 4, 6])
        self.assertEqual([v['next_output_ordinal'] for v in omissions], [1, 2, 2, 3])
        self.assertTrue(all(v['source_row_sha256'] == hashlib.sha256(missing).hexdigest() for v in omissions))
        self.assertEqual(ledger['sha256'], M.digest(path))
        self.assertEqual(ledger['records'], 4)
        self.assertEqual([(v['source_ordinal_start'], v['source_ordinal_end']) for v in proof['files']], [(1, 6), (7, 7)])
        self.assertEqual((self.root/'preserved/contig-000.vep.vcf').read_bytes(),
                         header + annotated(ROWS[0], 1) + annotated(ROWS[1], 4))
        self.assertEqual(self.serial.read_bytes(), original)
        self.assertTrue(proof['source_ordinal_occurrences_accounted'])
        self.assertFalse(proof['source_ordinal_occurrences_preserved'])
        self.assertFalse(proof['canonical_publication'])

    def test_a_completed_contig_with_only_recorded_nonvariants_is_header_only(self):
        missing = b'chr1\t10\t.\tT\t.\t.\tPASS\t.\tGT\t0/0\n'
        header = self.nonvariant_fixture([missing, missing, ROWS[2], ROWS[3]], [2, 3])
        proof = self.preserve()
        self.assertEqual(proof['records_by_contig'], {'chr1': 0, 'chr2': 1})
        self.assertEqual(proof['omitted_records_by_contig'], {'chr1': 2, 'chr2': 0})
        self.assertEqual((self.root/'preserved/contig-000.vep.vcf').read_bytes(), header)

    def test_unknown_missing_variant_still_fails_and_retains_omission_history(self):
        missing = b'chr1\t10\t.\tT\t.\t.\tPASS\t.\tGT\t0/0\n'
        self.nonvariant_fixture([missing, ROWS[0], ROWS[1], ROWS[2], ROWS[3]], [2, 3, 4])
        with self.assertRaisesRegex(ValueError, 'changed prefix occurrence'):
            self.preserve()
        self.assertFalse((self.root/'preserved/preserved-prefix.json').exists())
        ledger = self.root/'preserved/omitted-source-occurrences.jsonl'
        self.assertEqual([json.loads(line)['source_ordinal'] for line in ledger.read_text().splitlines()], [1])
        self.assertTrue((self.root/'preserved/contig-000.vep.vcf').exists())

    def test_nonvariant_omission_requires_pinned_source_and_effective_header(self):
        missing = b'chr1\t10\t.\tT\t.\t.\tPASS\t.\tGT\t0/0\n'
        self.nonvariant_fixture([missing, ROWS[0], ROWS[2], ROWS[3]], [1, 2, 3])
        self.vep_source_root = None
        with self.assertRaisesRegex(ValueError, 'changed prefix occurrence'):
            self.preserve()
        self.assertFalse((self.root/'preserved/preserved-prefix.json').exists())

    def test_nonvariant_source_pin_mismatch_has_no_ready_proof(self):
        self.nonvariant_fixture(ROWS, range(len(ROWS)))
        next(self.vep_source_root.rglob('Config.pm')).write_text('changed source\n')
        with self.assertRaisesRegex(ValueError, 'Pinned VEP116 omission source differs'):
            self.preserve()
        self.assertFalse((self.root/'preserved/preserved-prefix.json').exists())

    def test_allow_nonvariant_or_individual_flags_do_not_authorize_omission(self):
        for command in [b'vep --allow_non_variant', b'vep --allow_non_variant 1',
                        b'vep --individual SAMPLE', b'vep --individual_zyg all', b'vep --individual 0']:
            with self.subTest(command=command):
                header = b'##VEP="v116"\n##VEP-command-line=\'' + command + b"'\n"
                with mock.patch.object(M, 'VEP116_OMISSION_SOURCES', {}):
                    policy = M.vep116_omission_policy(header, self.root)
                self.assertTrue(policy['allow_non_variant'])
                self.assertFalse(policy['nonvariant_omissions_authorized'])

    def test_reference_and_dont_skip_flags_do_not_change_early_omission(self):
        for command in [b'vep --offline --use_transcript_ref --dont_skip',
                        b'vep --allow_non_variant 0 --check_ref --lookup_ref']:
            with self.subTest(command=command), mock.patch.object(M, 'VEP116_OMISSION_SOURCES', {}):
                policy = M.vep116_omission_policy(b'##VEP="v116.3"\n##VEP-command-line=\'' + command + b"'\n", self.root)
                self.assertFalse(policy['allow_non_variant'])
                self.assertTrue(policy['nonvariant_omissions_authorized'])

    def test_missing_duplicate_wrong_version_and_ambiguous_headers_fail_closed(self):
        valid = b"##VEP=\"v116\"\n##VEP-command-line='vep --offline'\n"
        for header in [b'', valid.replace(b'v116', b'v115'), valid + valid,
                       valid.replace(b'--offline', b'--allow_non_variant=0'),
                       valid.replace(b'--offline', b'--allow_non_variant 0 --allow_non_variant 0')]:
            with self.subTest(header=header), mock.patch.object(M, 'VEP116_OMISSION_SOURCES', {}):
                with self.assertRaises(ValueError):
                    M.vep116_omission_policy(header, self.root)

    def test_exact_parser_predicate_respects_structural_dispatch_and_perl_info_truth(self):
        cases = [(b'.', b'.', True), (b'.,A', b'.', True),
                 (b'.,<DEL>', b'.', False), (b'.,A.', b'.', False),
                 (b'.,A]chr2:2]', b'.', False), (b'.', b'SVTYPE=DEL', False),
                 (b'.', b'SVTYPE', True), (b'.', b'SVTYPE=0', True),
                 (b'.', b'SVTYPE=00', False), (b'.', b'SVTYPE==DEL', True),
                 (b'.', b'SVTYPE=DEL;SVTYPE=0', True), (b'.', b'SVTYPE=0;SVTYPE=DEL', False),
                 (b'T', b'.', False), (b'*', b'.', False), (b'<DEL>', b'.', False)]
        for alt, info, omitted in cases:
            fields = [b'chr1', b'10', b'id', b'T', alt, b'.', b'PASS', info]
            with self.subTest(alt=alt, info=info):
                self.assertEqual(M.vep116_nonvariant_omitted(fields, allow_non_variant=False), omitted)
                self.assertFalse(M.vep116_nonvariant_omitted(fields, allow_non_variant=True))

    def test_corrupted_omission_ledger_has_no_ready_proof(self):
        missing = b'chr1\t10\t.\tT\t.\t.\tPASS\t.\tGT\t0/0\n'
        self.nonvariant_fixture([missing, ROWS[0], ROWS[2], ROWS[3]], [1, 2, 3])
        original_digest = M.digest
        def corrupt_ledger(path):
            if Path(path).name == 'omitted-source-occurrences.jsonl':
                with Path(path).open('ab') as f:
                    f.write(b'corrupted\n')
            return original_digest(path)
        with mock.patch.object(M, 'digest', side_effect=corrupt_ledger), self.assertRaisesRegex(ValueError, 'ledger checksum differs'):
            self.preserve()
        self.assertFalse((self.root/'preserved/preserved-prefix.json').exists())


if __name__ == '__main__':
    unittest.main()
