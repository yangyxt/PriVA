#!/usr/bin/env python3
"""Synthetic selection-view mechanics, never biological qualification."""
import gzip
import copy
import json
from pathlib import Path
import sys
import tempfile
import unittest

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / 'bin'))
import cadd_select_view as V


class SelectionViewTests(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory()
        self.root = Path(self.temporary.name)
        self.source, self.output = self.root / 'source.vcf', self.root / 'view.vcf.gz'
        self.columns = V.S.REQUIRED_CSQ + ['CADD_PICK_ORDINAL', 'CADD_PICK_COUNT', 'CADD_PICK_CONTEXT']
        self.context = 'a' * 64

    def tearDown(self):
        self.temporary.cleanup()

    def write(self, orders=(2, 0, 1), counts=(2, 2, 2), context=None):
        annotations = []
        for index, (order, count) in enumerate(zip(orders, counts)):
            row = dict.fromkeys(self.columns, '')
            row.update(Allele='G', Consequence='missense_variant', Feature_type='Transcript',
                       Gene='GENE', Feature='TRANSCRIPT' + str(index), CADD_PICK_ORDINAL=str(order),
                       CADD_PICK_COUNT=str(count), CADD_PICK_CONTEXT=context or self.context)
            annotations.append('|'.join(row[c] for c in self.columns))
        self.source.write_text('##fileformat=VCFv4.2\n##INFO=<ID=CSQ,Number=.,Type=String,Description="Format: '
                               + '|'.join(self.columns) + '">\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n'
                               + '1\t100\t.\tA\tG\t.\tPASS\tKEEP=yes;CSQ=' + ','.join(annotations) + '\n')

    def test_native_order_and_original_source_retained(self):
        self.write()
        original = self.source.read_bytes()
        self.assertEqual(V.select_rows(self.source, self.output, self.context), 1)
        self.assertEqual(self.source.read_bytes(), original)
        self.assertEqual([r[-1] for r in V.S.expected_features(self.output)], ['TRANSCRIPT2', 'TRANSCRIPT0'])
        self.assertIn('KEEP=yes;', gzip.open(self.output, 'rt').read())
        self.assertEqual(list(V.S.variant_keys(self.source)), list(V.S.variant_keys(self.output)))

    def test_incomplete_duplicate_and_missing_plugin_rows_rejected(self):
        for orders, counts, message in [((1, 0), (2, 2), 'incomplete'),
                                        ((1, 1), (2, 2), 'Duplicate'),
                                        (('', 0), (2, 2), 'invalid'),
                                        ((1, 2), (2, 3), 'incomplete')]:
            with self.subTest(orders=orders, counts=counts):
                self.write(orders, counts)
                with self.assertRaisesRegex(V.S.SharedVepError, message):
                    V.select_rows(self.source, self.output, self.context)

    def test_wrong_native_context_rejected(self):
        self.write(context='b' * 64)
        with self.assertRaisesRegex(V.S.SharedVepError, 'context differs'):
            V.select_rows(self.source, self.output, self.context)

    def profile(self):
        context = {'adapter_id': V.ADAPTER, 'vep_release': '113', 'perl_version': 'v5.fixture',
                   'plugin_sha256': '1' * 64, 'output_factory_sha256': '2' * 64,
                   'perl_hash_seed': '0', 'perl_perturb_keys': '0', 'pick_order': ['canonical', 'length'],
                   'cadd_transcript_sources': ['Ensembl'],
                   'annotation_options': dict(V.S.NATIVE_REQUIRED_OPTIONS, assembly='GRCh37')}
        policy = dict(copy.deepcopy(V.S.SHARED_TRANSCRIPT_POLICY), cadd_transcript_sources=['Ensembl'])
        return {'profile_kind': 'adapted', 'profile_id': 'synthetic-native-view-only',
                'contract': {'assembly': 'GRCh37', 'resource_manifest_sha256': '0' * 64},
                'annotation': {'shared_vep': {'adapter_id': V.ADAPTER, 'vep_release': '113',
                    'transcript_release': '113', 'esm_transcript_release': '113',
                    'transcript_policy': policy,
                    'selection_context': context}}}

    def test_fixed_shared_profile_mismatches_rejected(self):
        profile = self.profile()
        self.assertEqual(len(V.S.native_context(profile)), 64)
        changes = [('profile_kind', 'reference'), ('merged', ''), ('use_transcript_ref', ''),
                   ('assembly', 'GRCh38'), ('sift', ''), ('regulatory', '')]
        for key, value in changes:
            bad = copy.deepcopy(profile)
            if key == 'profile_kind':
                bad[key] = value
            else:
                bad['annotation']['shared_vep']['selection_context']['annotation_options'][key] = value
            with self.subTest(key=key), self.assertRaises(V.S.SharedVepError):
                V.S.native_context(bad)
        bad = copy.deepcopy(profile)
        bad['annotation']['shared_vep']['transcript_policy']['cadd_transcript_sources'] = ['Ensembl', 'RefSeq']
        with self.assertRaises(V.S.SharedVepError):
            V.S.native_context(bad)

    def test_completed_view_records_qualified_context_and_coverage(self):
        profile = self.profile()
        settings = profile['annotation']['shared_vep']
        self.context = V.S.native_context(profile)
        self.write()
        proof = self.root / 'synthetic-proof.json'
        proof.write_text('{"scope":"synthetic-wiring-only"}')
        qualification = self.root / 'qualification.json'
        qualification.write_text(json.dumps({'schema_version': 1, 'status': 'qualified', 'assembly': 'GRCh37',
            'vep_release': '113', 'transcript_release': '113', 'selection': 'cadd-per-gene',
            'adapter_id': V.ADAPTER, 'stub': True, 'selection_context_sha256': self.context,
            'transcript_policy': settings['transcript_policy'],
            'checks': {'per_gene_selection': 'exact', 'feature_id_joins': 'exact',
                       'variant_representation': 'exact', 'missingness': 'preserved'},
            'evidence': {'path': str(proof), 'sha256': V.S.sha256(proof)}}))
        settings['qualification'] = {'path': str(qualification), 'sha256': V.S.sha256(qualification)}
        bundle = self.root / 'bundle'
        V.build(self.source, profile, bundle, stub=True)
        receipt = V.S.json_file(bundle / 'manifest.json')
        self.assertEqual(receipt['selection_context_sha256'], self.context)
        self.assertEqual(receipt['expected_annotation_rows'], 2)
        # Even an internally hash-consistent receipt must not drop variants.
        original = self.source.read_text()
        self.source.write_text(original + original.splitlines()[-1] + '\n')
        receipt['source_vep']['sha256'] = V.S.sha256(self.source)
        (bundle / 'manifest.json').write_text(json.dumps(receipt))
        with self.assertRaisesRegex(V.S.SharedVepError, 'coverage'):
            V.S.verify(bundle, profile, stub=True)


if __name__ == '__main__':
    unittest.main()
