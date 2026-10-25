#!/usr/bin/env python3
"""Small synthetic legacy-versus-checkpoint tests; no real models or caches.

The caller must supply a fresh, explicitly authorized filesystem fixture root.
All fixtures, interrupted attempts and result metadata remain there for review.
These tests do not establish real-model, installation or family qualification.
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
import sys
import time
import unittest
from unittest.mock import patch

sys.dont_write_bytecode = True
SOURCE = Path(__file__).absolute().parents[1] / 'scripts/alphamissense_domain_cache.py'
spec = importlib.util.spec_from_file_location('domain_checkpoint_under_test', SOURCE)
m = importlib.util.module_from_spec(spec); spec.loader.exec_module(m)
FIXTURE_ROOT = None


def annotation(gene, *, domains='PANTHER:P1:SF1&PANTHER:P1', transcript='ENST0001', exon='1/3', position='10/100'):
    return '|'.join(('C', gene, 'Transcript', transcript, exon, position, domains, ''))


class DomainCheckpoints(unittest.TestCase):
    def setUp(self):
        self.folder = FIXTURE_ROOT / self._testMethodName
        self.folder.mkdir()
        self.cache = self.folder / 'cache'; self.cache.mkdir()
        self.vcf = self.cache / 'input.vcf.gz'
        self.output = self.cache / 'domain.pkl'
        self.module = m.producer()
        self.contigs = ['chr2', 'chrX', 'chr1', 'empty', 'chrY']
        self.make_vcf()

    def make_vcf(self, *, malformed=False, missing=False):
        plain = self.cache / 'input.vcf'
        header = ['##fileformat=VCFv4.2']
        header += ['##contig=<ID=' + chrom + ',length=1000>' for chrom in self.contigs]
        header += [
            '##INFO=<ID=AM_PATHOGENICITY,Number=1,Type=Float,Description="Synthetic score">',
            '##INFO=<ID=AF_grpmax_joint,Number=A,Type=Float,Description="Synthetic AF">',
            '##INFO=<ID=PVAR,Number=1,Type=String,Description="Synthetic PVAR">',
            '##INFO=<ID=CSQ,Number=.,Type=String,Description="Synthetic. Format: Allele|Gene|Feature_type|Feature|EXON|Protein_position|DOMAINS|DUMMY">',
            '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO',
        ]
        records = [
            ('chr2', 10, 'p.R10C', 0.2, 0.1, ','.join([annotation('GENEA')] * 2)),
            ('chrX', 11, 'p.R11C', 0.4, 0.2, annotation('PAR', transcript='ENSTPAR', position='11/100')),
            ('chr1', 12, 'p.R10C', 0.9, 0.8, annotation('GENEA', exon='2/3', position='12/100')),
            ('chr1', 13, 'p.R12G', 0.6, None, annotation('OTHER', domains='Pfam:PF0001', position='12/100')),
            ('chrY', 14, 'p.R11G', 0.3, 0.9, annotation('PAR', transcript='ENSTPAR', exon='3/3', position='14/100')),
        ]
        body = []
        for chrom, pos, pvar, score, af, csq in records:
            info = ['AM_PATHOGENICITY=' + str(score), 'PVAR=' + pvar, 'CSQ=' + csq]
            if af is not None: info.append('AF_grpmax_joint=' + str(af))
            if missing and chrom == 'chr2': info = [value for value in info if not value.startswith('PVAR=')]
            if malformed and chrom == 'chr2': info = ['AM_PATHOGENICITY=0.2', 'PVAR=p.R10C', 'CSQ=C']
            body.append('\t'.join((chrom, str(pos), '.', 'A', 'C', '.', 'PASS', ';'.join(info))))
        plain.write_text('\n'.join(header + body) + '\n')
        self.module.pysam.tabix_compress(str(plain), str(self.vcf), force=True)
        self.module.pysam.tabix_index(str(self.vcf), preset='vcf', force=True)

    def plan(self):
        return m.create_plan(input_vcf=self.vcf, output_pickle=self.output, cache_root=self.cache)

    def collect_all(self, plan_path):
        plan = m.load_plan(plan_path)
        return [m.collect_chromosome(plan_path=plan_path, ordinal=i) for i in range(len(plan['recipe']['contigs']))]

    def legacy(self):
        result = {}
        # Exact original parallel reduction: all raw contigs, header order,
        # unchanged merge/finalize/to_plain functions. No serial-mode oracle.
        parts = [self.module._collect_worker((str(self.vcf), chrom)) for chrom in self.contigs]
        for part in parts:
            self.module._merge_scores(result, part)
        self.module._finalize_scores_dict(result)
        return self.module.to_plain_dict(result)

    def loaded(self):
        with self.output.open('rb') as stream:
            return pickle.load(stream)

    def test_full_original_reduction_equals_checkpointed_ordered_gene_merge(self):
        expected = self.legacy()
        result = m.build(input_vcf=self.vcf, output_pickle=self.output, cache_root=self.cache)
        self.assertFalse(result['reused'])
        actual = self.loaded(); m.assert_equal(expected, actual)
        self.assertEqual(list(actual), ['GENEA', 'PAR', 'OTHER'])
        self.assertEqual(actual['GENEA']['PANTHER']['P1']['ENST0001']['exons'], {'1/3', '2/3'})
        self.assertEqual(actual['PAR']['PANTHER']['P1']['ENSTPAR']['aa_min'], 11)
        self.assertEqual(actual['PAR']['PANTHER']['P1']['ENSTPAR']['aa_max'], 14)

    def test_duplicate_csq_domains_first_wins_af_max_and_residue_multiplicity(self):
        import numpy as np
        m.build(input_vcf=self.vcf, output_pickle=self.output, cache_root=self.cache)
        actual = self.loaded()['GENEA']['PANTHER']['P1']
        legacy = self.legacy()['GENEA']['PANTHER']['P1']; m.assert_equal(legacy, actual)
        self.assertEqual(len(actual['distribution']), 1)
        # Two duplicate CSQs times three parent-domain visits, then chr1's
        # three visits. Nothing may silently deduplicate that original list.
        self.assertEqual(len(actual['gnomAD_residue_distribution']['p.R10']), 9)
        self.assertLess(actual['distribution'][0], actual['max_distribution'][0])
        self.assertEqual(actual['distribution'].tolist(), [float(np.float32(0.2))])
        self.assertEqual(actual['min_distribution'].tolist(), [float(np.float32(0.2))])
        self.assertEqual(actual['max_distribution'].tolist(), [float(np.float32(0.9))])
        self.assertEqual(actual['gnomAD_distribution'].tolist(), [float(np.float32(0.8))])
        self.assertEqual(actual['gnomAD_residue_distribution']['p.R10'],
                         [float(np.float32(0.1))] * 6 + [float(np.float32(0.8))] * 3)

    def test_literal_header_order_empty_chromosome_and_safe_ordinal_names(self):
        path = self.plan(); plan = m.load_plan(path)
        self.assertEqual(plan['recipe']['contigs'], self.contigs)
        value = m.collect_chromosome(plan_path=path, ordinal=3)
        receipt = m.read_json(value['receipt'])
        self.assertEqual(m.read_json(receipt['genes']['path']), [])
        self.assertEqual(Path(receipt['raw']['path']).parent.name, 'chromosome-00003')
        with Path(receipt['raw']['path']).open('rb') as stream: self.assertEqual(pickle.load(stream), {})

    def test_content_identical_source_paths_do_not_change_generation(self):
        path = self.plan(); original = m.source_state
        def elsewhere():
            value = original()
            for item in value.values(): item['path'] = '/different/frozen/snapshot/' + Path(item['path']).name
            return value
        with patch.object(m, 'source_state', side_effect=elsewhere):
            self.assertEqual(self.plan(), path)

    def test_completed_checkpoint_reuse_retains_exact_raw_bytes_and_identity(self):
        path = self.plan(); first = m.collect_chromosome(plan_path=path, ordinal=0)
        raw = m.read_json(first['receipt'])['raw']; before = Path(raw['path']).read_bytes()
        second = m.collect_chromosome(plan_path=path, ordinal=0)
        self.assertTrue(second['reused']); self.assertEqual(Path(raw['path']).read_bytes(), before)
        self.assertEqual(m.identity(raw['path']), raw['identity'])

    def test_corrupt_owned_checkpoint_is_retained_and_recollected(self):
        path = self.plan(); first = m.collect_chromosome(plan_path=path, ordinal=0)
        receipt = m.read_json(first['receipt']); raw = Path(receipt['raw']['path'])
        with raw.open('ab') as stream: stream.write(b'corrupt')
        corrupt = raw.read_bytes()
        second = m.collect_chromosome(plan_path=path, ordinal=0)
        self.assertFalse(second['reused']); self.assertEqual(raw.read_bytes(), corrupt)
        self.assertNotEqual(m.read_json(second['receipt'])['raw']['path'], str(raw))
        self.assertEqual(len(list(raw.parent.glob('complete.invalid-*.json'))), 1)

    def test_unfinished_attempt_is_preserved_and_not_reused(self):
        path = self.plan(); folder, ready = m.chromosome_paths(m.load_plan(path), 0)
        folder.mkdir(); partial = folder / 'old.raw.pkl.partial-retained'; partial.write_bytes(b'unfinished')
        result = m.collect_chromosome(plan_path=path, ordinal=0)
        self.assertFalse(result['reused']); self.assertEqual(partial.read_bytes(), b'unfinished')
        self.assertTrue(ready.exists())

    def test_missing_checkpoint_cannot_publish_any_output(self):
        path = self.plan(); m.collect_chromosome(plan_path=path, ordinal=0)
        with self.assertRaises(m.InvalidCache): m.merge(plan_path=path)
        self.assertFalse(self.output.exists())

    def test_malformed_original_indexerror_is_not_silently_skipped(self):
        self.make_vcf(malformed=True); path = self.plan()
        with self.assertRaises(IndexError): self.module._collect_worker((str(self.vcf), 'chr2'))
        with self.assertRaises(IndexError): m.collect_chromosome(plan_path=path, ordinal=0)
        self.assertFalse(m.chromosome_paths(m.load_plan(path), 0)[1].exists())

    def test_missing_original_pvar_behavior_is_preserved(self):
        self.make_vcf(missing=True); expected = self.legacy()
        m.build(input_vcf=self.vcf, output_pickle=self.output, cache_root=self.cache)
        m.assert_equal(expected, self.loaded())

    def test_input_change_during_collection_prevents_completion(self):
        path = self.plan(); original = self.module.DomainAMScoreCollector.collect_scores
        def changed(collector, region=None):
            original(collector, region=region)
            os.utime(self.vcf, ns=(self.vcf.stat().st_atime_ns, self.vcf.stat().st_mtime_ns + 1))
        with patch.object(m, 'producer', return_value=self.module), patch.object(self.module.DomainAMScoreCollector, 'collect_scores', changed):
            with self.assertRaisesRegex(ValueError, 'changed during'): m.collect_chromosome(plan_path=path, ordinal=0)
        plan = m.read_json(path)
        self.assertFalse(m.chromosome_paths(plan, 0)[1].exists())

    def test_source_identity_change_during_collection_prevents_completion(self):
        path = self.plan(); original = m.source_state; counter = [0]
        def observed():
            counter[0] += 1; value = original()
            if counter[0] >= 3: value['producer']['identity']['st_mtime_ns'] = '1'
            return value
        with patch.object(m, 'source_state', side_effect=observed):
            with self.assertRaisesRegex(ValueError, 'changed during'): m.collect_chromosome(plan_path=path, ordinal=0)
        self.assertFalse(m.chromosome_paths(m.load_plan(path), 0)[1].exists())

    def test_unknown_canonical_is_preserved_before_checkpoint_directory_writes(self):
        self.output.write_bytes(b'user supplied data')
        with self.assertRaisesRegex(ValueError, 'unowned'): self.plan()
        self.assertEqual(self.output.read_bytes(), b'user supplied data')
        self.assertFalse(Path(str(self.output) + '.chromosomes').exists())

    def test_input_and_index_aliases_fail_before_writes(self):
        os.link(str(self.vcf) + '.tbi', self.output)
        before = self.output.read_bytes()
        with self.assertRaisesRegex(ValueError, 'aliases'): self.plan()
        self.assertEqual(self.output.read_bytes(), before)
        self.assertFalse(Path(str(self.output) + '.chromosomes').exists())

    def test_symlink_destination_escape_is_preserved(self):
        other = self.folder / 'foreign'; other.mkdir()
        link = self.cache / 'linked'; link.symlink_to(other, target_is_directory=True)
        with self.assertRaisesRegex(ValueError, 'aliases outside|linked destination'):
            m.create_plan(input_vcf=self.vcf, output_pickle=link / 'output.pkl', cache_root=self.cache)
        self.assertEqual(list(other.iterdir()), [])

    def test_ready_receipt_failure_leaves_pending_and_retry_retains_exact_payload(self):
        path = self.plan(); self.collect_all(path)
        original = m.atomic_json; ready = Path(str(self.output) + '.complete.json')
        def interrupted(target, value):
            if Path(target) == ready: raise OSError('synthetic receipt publication failure')
            return original(target, value)
        with patch.object(m, 'atomic_json', side_effect=interrupted):
            with self.assertRaisesRegex(OSError, 'synthetic'): m.merge(plan_path=path)
        before = self.output.read_bytes()
        self.assertTrue(Path(str(self.output) + '.publish.pending.json').exists())
        self.assertIsNone(m.output_ready(m.load_plan(path)))
        m.merge(plan_path=path)
        m.assert_equal(self.legacy(), self.loaded())
        retained = list(Path(str(self.output) + '.chromosomes').glob('retained-*/domain.pkl'))
        self.assertEqual([item.read_bytes() for item in retained], [before])
        self.assertFalse(Path(str(self.output) + '.publish.pending.json').exists())

    def test_checkpoint_write_failure_retains_partial_and_no_ready(self):
        path = self.plan(); original = m.pickle.dump
        def interrupted(value, stream, **kwargs):
            stream.write(b'partial'); raise OSError('synthetic pickle failure')
        with patch.object(m.pickle, 'dump', side_effect=interrupted):
            with self.assertRaisesRegex(OSError, 'synthetic'): m.collect_chromosome(plan_path=path, ordinal=0)
        folder, ready = m.chromosome_paths(m.load_plan(path), 0)
        self.assertFalse(ready.exists()); self.assertEqual(len(list(folder.glob('*.partial-*'))), 1)
        m.collect_chromosome(plan_path=path, ordinal=0)
        self.assertTrue(ready.exists()); self.assertEqual(len(list(folder.glob('*.partial-*'))), 1)

    def test_pickle_protocol_five_buffer_roundtrip_and_complete_hash(self):
        import numpy as np
        value = dict(scores=np.arange(32768, dtype=np.float64))
        observed = []
        original = m.HashWriter
        class ObservedWriter(original):
            def write(self, buffer):
                observed.append(isinstance(buffer, pickle.PickleBuffer))
                return super().write(buffer)
        target = self.cache / 'large-array.pkl'
        with patch.object(m, 'HashWriter', ObservedWriter):
            result = m.write_pickle(target, value)
        self.assertTrue(any(observed), 'Regression must execute the actual PickleBuffer path')
        self.assertEqual(result['sha256'], hashlib.sha256(target.read_bytes()).hexdigest())
        with target.open('rb') as stream:
            m.assert_equal(value, pickle.load(stream))
        self.assertEqual(list(self.cache.glob('large-array.pkl.partial-*')), [])

    def test_typed_buffer_write_counts_bytes_instead_of_elements(self):
        from array import array
        from io import BytesIO
        value = array('I', [0, 1, 256, 65536])
        with memoryview(value) as view:
            self.assertNotEqual(len(view), view.nbytes)
            stream = BytesIO()
            writer = m.HashWriter(stream)
            self.assertEqual(writer.write(view), view.nbytes)
        self.assertEqual(stream.getvalue(), value.tobytes())
        self.assertEqual(writer.hash.hexdigest(), hashlib.sha256(value.tobytes()).hexdigest())

    def test_actual_pickle_buffer_short_write_or_io_failure_never_publishes(self):
        import numpy as np
        value = np.arange(32768, dtype=np.float64)
        original = m.HashWriter
        for mode in ('short', 'io-error'):
            target = self.cache / (mode + '.pkl')
            observed = []
            class FaultyStream:
                def __init__(self, stream):
                    self.stream = stream
                def write(self, buffer):
                    if isinstance(buffer, pickle.PickleBuffer):
                        observed.append(True)
                        if mode == 'io-error':
                            raise OSError('synthetic buffer write failure')
                        with memoryview(buffer) as view:
                            return self.stream.write(view.cast('B')[:-1])
                    return self.stream.write(buffer)
            class FaultyWriter(original):
                def __init__(self, stream):
                    super().__init__(FaultyStream(stream))
            error, message = ((ValueError, 'Short checkpoint write') if mode == 'short'
                              else (OSError, 'synthetic buffer write failure'))
            with self.subTest(mode=mode), patch.object(m, 'HashWriter', FaultyWriter):
                with self.assertRaisesRegex(error, message):
                    m.write_pickle(target, value)
            self.assertTrue(observed, 'Failure must occur in the actual PickleBuffer path')
            self.assertFalse(target.exists())
            retained = list(self.cache.glob(mode + '.pkl.partial-*'))
            self.assertEqual(len(retained), 1)
            before = retained[0].read_bytes()
            m.write_pickle(target, value)
            self.assertEqual(retained[0].read_bytes(), before)
            with target.open('rb') as stream:
                m.assert_equal(value, pickle.load(stream))

    def test_full_model_reuse_and_explicit_worker_allocation_cap(self):
        m.build(input_vcf=self.vcf, output_pickle=self.output, cache_root=self.cache)
        before = m.identity(self.output)
        self.assertTrue(m.build(input_vcf=self.vcf, output_pickle=self.output, cache_root=self.cache)['reused'])
        self.assertEqual(m.identity(self.output), before)
        with patch.dict(os.environ, SLURM_CPUS_PER_TASK='1'):
            with self.assertRaisesRegex(ValueError, 'CPU allocation'):
                m.build(input_vcf=self.vcf, output_pickle=self.output, cache_root=self.cache, workers=2)

    def test_reuse_binds_every_ordered_completion_receipt_without_raw_rehash(self):
        m.build(input_vcf=self.vcf, output_pickle=self.output, cache_root=self.cache)
        plan = m.load_plan(self.plan())
        ready = Path(str(self.output) + '.complete.json')
        original = m.read_json(ready)
        cp_path = Path(original['checkpoint_receipts'][0]['path'])
        cp_original = m.read_json(cp_path)
        original_digest = m.digest
        def no_raw_rehash(path):
            self.assertFalse(str(path).endswith('.raw.pkl'), 'Reuse must not rehash raw chromosome payloads')
            return original_digest(path)
        with patch.object(m, 'digest', side_effect=no_raw_rehash):
            self.assertIsNotNone(m.output_ready(plan))
            for mode in ('missing', 'empty', 'order', 'sha', 'path', 'generation', 'ordinal', 'chromosome', 'complete'):
                value = copy.deepcopy(original)
                if mode == 'missing': value.pop('checkpoint_receipts')
                elif mode == 'empty': value['checkpoint_receipts'] = []
                elif mode == 'order': value['checkpoint_receipts'].reverse()
                elif mode == 'sha': value['checkpoint_receipts'][0]['sha256'] = '0' * 64
                elif mode == 'path': value['checkpoint_receipts'][0]['path'] = value['checkpoint_receipts'][1]['path']
                else:
                    changed = copy.deepcopy(cp_original)
                    changed[mode] = {'generation': 'wrong', 'ordinal': 1, 'chromosome': 'wrong', 'complete': False}[mode]
                    m.atomic_json(cp_path, changed)
                    value['checkpoint_receipts'][0]['sha256'] = m.digest(cp_path)
                m.atomic_json(ready, value)
                with self.subTest(mode=mode): self.assertIsNone(m.output_ready(plan))
                m.atomic_json(cp_path, cp_original)
            m.atomic_json(ready, original)
            self.assertIsNotNone(m.output_ready(plan))

    def test_exact_comparator_checks_nan_signed_zero_dtype_and_order(self):
        import numpy as np
        values = dict(a=np.array([0.0, -0.0, float('nan')], dtype=np.float64),
                      b=[float('nan'), -0.0, float('inf')], exons={'1/3', '2/3'})
        m.assert_equal(values, copy.deepcopy(values))
        for change in ('dtype', 'zero', 'order', 'set', 'list'):
            other = copy.deepcopy(values)
            if change == 'dtype': other['a'] = other['a'].astype(np.float32)
            elif change == 'zero': other['a'][1] = 0.0
            elif change == 'order': other = dict(reversed(list(other.items())))
            elif change == 'set': other['exons'].add('3/3')
            else: other['b'].append(0.0)
            with self.subTest(change=change), self.assertRaises(ValueError): m.assert_equal(values, other)

    def test_later_full_chromosome_comparison_entrypoint_against_independent_legacy_collection(self):
        path = self.plan(); m.collect_chromosome(plan_path=path, ordinal=0)
        result = m.compare_chromosome(plan_path=path, chromosome='chr2')
        self.assertTrue(result['equal']); self.assertFalse(result['scientific_qualification'])


def main(argv=None):
    global FIXTURE_ROOT
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest='procedure', required=True)
    child = sub.add_parser('test'); child.add_argument('--fixture-root', type=Path, required=True)
    args = parser.parse_args(argv)
    FIXTURE_ROOT = args.fixture_root.absolute()
    require_root = Path('/stornext/snfs170/next-gen/scratch/yangyxt_dt/projects/PriVA_upgrade/tmp/alphamissense-domain-cache-fixtures-20261024-v1')
    m.require(FIXTURE_ROOT.is_relative_to(require_root) and FIXTURE_ROOT != require_root and not FIXTURE_ROOT.exists(),
              'Use a fresh child of the exact authorized fixture prefix')
    FIXTURE_ROOT.mkdir(parents=True)
    before = m.source_state()
    test_sha = m.digest(Path(__file__))
    started = time.monotonic()
    print(json.dumps(dict(event='fixture_start', native_step=os.environ.get('SLURM_JOB_ID', '') + '.' + os.environ.get('SLURM_STEP_ID', ''),
                          source=before, tests_sha256=test_sha, fixture_root=str(FIXTURE_ROOT))), flush=True)
    result = unittest.TextTestRunner(verbosity=2, stream=sys.stdout).run(unittest.defaultTestLoader.loadTestsFromTestCase(DomainCheckpoints))
    unchanged = m.source_state() == before and m.digest(Path(__file__)) == test_sha
    value = dict(tests=result.testsRun, failures=len(result.failures), errors=len(result.errors), skipped=len(result.skipped),
                 source_unchanged=unchanged, source=before, tests_sha256=test_sha, wall_seconds=time.monotonic() - started,
                 native_step=os.environ.get('SLURM_JOB_ID', '') + '.' + os.environ.get('SLURM_STEP_ID', ''),
                 finished_utc=datetime.now(timezone.utc).isoformat(), installation_accepted=False,
                 scientific_qualification=False, real_model_validation=False, family_qualification=False)
    m.atomic_json(FIXTURE_ROOT / 'result.json', value)
    print(json.dumps(value), flush=True)
    return 0 if result.wasSuccessful() and unchanged else 1


if __name__ == '__main__':
    raise SystemExit(main())
