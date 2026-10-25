"""Real indexed tiny inputs: chromosome failures must abort parent publication."""
import ast
import copy
import importlib
import multiprocessing as mp
import os
from pathlib import Path
import pickle
import subprocess
import sys
import tempfile
import unittest
from unittest.mock import patch

import pysam

SCRIPTS = Path(__file__).resolve().parents[1] / 'scripts'
sys.path.insert(0, str(SCRIPTS))
NAMES = ('am_pick_missense_constrained_segments', 'stat_aachange_clinvar', 'stat_exon_patho_afs')
MODULES = [importlib.import_module(name) for name in NAMES]
FIELDS = 'Gene|Feature|Consequence|Protein_position|HGVSp|Feature_type|EXON|INTRON|HGVSc|Unused'
REAL_VARIANT = pysam.VariantFile


class FaultFile:
    """Yield one real record, then inject an unexpected chromosome read error."""
    def __init__(self, path, error):
        self.file = REAL_VARIANT(path)
        self.error = error

    def fetch(self, chromosome):
        yield next(self.file.fetch(chromosome))
        raise self.error('injected chromosome read failure')

    def close(self):
        self.file.close()


class WorkerFailures(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)
        raw = self.root / 'tiny.vcf'
        header = '##fileformat=VCFv4.2\n##contig=<ID=1>\n##contig=<ID=2>\n'
        for key, number, kind in [('AM_PATHOGENICITY','1','String'), ('TRANSCRIPT','1','String'),
                                  ('PVAR','1','String'), ('CLNSIG','.','String'),
                                  ('CLNREVSTAT','.','String'), ('AF_grpmax_joint','.','String')]:
            header += f'##INFO=<ID={key},Number={number},Type={kind},Description="fixture">\n'
        header += f'##INFO=<ID=CSQ,Number=.,Type=String,Description="Format: {FIELDS}">\n'
        header += '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n'
        rows = []
        for chrom in ('1','2'):
            for pos in (10,20):
                info = (f'AM_PATHOGENICITY=0.8;TRANSCRIPT=ENST{chrom}.1;PVAR=A{pos}C;'
                        'CLNSIG=Pathogenic;CLNREVSTAT=reviewed_by_expert_panel;AF_grpmax_joint=0.02;'
                        f'CSQ=G{chrom}|ENST{chrom}|missense_variant|{pos}|p.A{pos}C|Transcript|1/2||c.1A>C|x')
                rows.append(f'{chrom}\t{pos}\t.\tA\tC\t.\tPASS\t{info}\n')
        raw.write_text(header + ''.join(rows))
        self.vcf = str(self.root / 'tiny.vcf.gz')
        pysam.tabix_compress(str(raw), self.vcf, force=True)
        pysam.tabix_index(self.vcf, preset='vcf', force=True)

    def worker(self, index, chromosome='1', legacy=False):
        module = MODULES[index]
        if index == 0:
            fn = module._parse_chromosome
            args = ((self.vcf, chromosome, 0, 1, 2),)
        elif index == 1:
            collector = module.AAClinvarCollector(self.vcf)
            fn = module.AAClinvarCollector.process_chromosome
            args = (collector, chromosome, collector.get_csq_indices())
        else:
            fn = module.process_chromosome
            args = (chromosome, self.vcf, {'Transcript':1,'EXON':6}, {'reviewed_by_expert_panel'})
        if legacy:
            # Reconstruct only the former outer exception handler. When the
            # frozen source is supplied, require exact AST equality to it.
            source = ast.parse(Path(module.__file__).read_text())
            owner = next(n for n in source.body if isinstance(n, ast.ClassDef)) if index == 1 else source
            node = next(n for n in owner.body if isinstance(n, ast.FunctionDef) and n.name == fn.__name__)
            tree = ast.Module(body=[copy.deepcopy(node)], type_ignores=[])
            outer = next(node for node in tree.body[0].body if isinstance(node, ast.Try))
            old = 'pass' if index == 0 else ('return ({}, {})' if index == 1 else 'return {}')
            if index == 0:
                outer.handlers[-1].body.pop()
            else:
                outer.handlers[-1].body[-1] = ast.parse(old).body[0]
            baseline = os.environ.get('PRIVA_WORKER_BASELINE')
            if baseline:
                frozen = ast.parse((Path(baseline) / 'scripts' / (NAMES[index] + '.py')).read_text())
                owner = next(n for n in frozen.body if isinstance(n, ast.ClassDef)) if index == 1 else frozen
                previous = next(n for n in owner.body if isinstance(n, ast.FunctionDef) and n.name == fn.__name__)
                self.assertEqual(ast.dump(tree.body[0]), ast.dump(previous))
            namespace = dict(module.__dict__)
            exec(compile(ast.fix_missing_locations(tree), '<verified former worker>', 'exec'), namespace)
            fn = namespace[fn.__name__]
        return fn(*args)

    def test_success_values_order_and_pickle_differential(self):
        for index in range(3):
            for chromosome in ('1','2'):
                with self.subTest(producer=NAMES[index], chromosome=chromosome):
                    actual = self.worker(index, chromosome)
                    self.assertTrue(actual)
                    self.assertEqual(pickle.dumps(actual), pickle.dumps(self.worker(index, chromosome, True)))

    def test_expected_record_keyerror_and_valueerror_still_skip(self):
        for index in range(3):
            for error in (KeyError, ValueError):
                class BadInfo(dict):
                    def __contains__(self, key): raise error('expected record skip')
                    def __getitem__(self, key): raise error('expected record skip')
                class BadRecord:
                    info = BadInfo()
                    chrom = '1'
                    pos = 9
                class SkipFile(FaultFile):
                    def fetch(self, chromosome):
                        yield BadRecord()
                        yield from self.file.fetch(chromosome)
                with self.subTest(producer=NAMES[index], error=error.__name__):
                    expected = self.worker(index)
                    with patch.object(pysam, 'VariantFile', lambda path: SkipFile(path, error)):
                        # AA header lookup needs the real handle, outside patch.
                        if index == 1:
                            with patch.object(MODULES[1].AAClinvarCollector, 'get_csq_indices', return_value=self.indices):
                                self.assertEqual(expected, self.worker(index))
                        else:
                            self.assertEqual(expected, self.worker(index))

    @property
    def indices(self):
        with patch.object(pysam, 'VariantFile', REAL_VARIANT):
            return MODULES[1].AAClinvarCollector(self.vcf).get_csq_indices()

    def test_unexpected_partial_worker_failure_is_not_a_result(self):
        for index in range(3):
            for error in (OSError, RuntimeError):
                with self.subTest(producer=NAMES[index], error=error.__name__):
                    indices = self.indices
                    with patch.object(pysam, 'VariantFile', lambda path: FaultFile(path, error)), \
                         patch.object(MODULES[1].AAClinvarCollector, 'get_csq_indices', return_value=indices):
                        self.worker(index, legacy=True)  # Original silently returns partial/empty.
                        with self.assertRaises(error): self.worker(index)

    def test_real_pool_child_failure_preserves_existing_outputs(self):
        for index in range(3):
            for error in (OSError, RuntimeError):
                with self.subTest(producer=NAMES[index], error=error.__name__):
                    out = self.root / 'output.pkl'
                    splice = self.root / 'splice.pkl'
                    out.write_bytes(b'previous output'); splice.write_bytes(b'previous splice')
                    before = [(p.read_bytes(), p.stat().st_mtime_ns) for p in (out, splice)]
                    def opened(path):
                        return REAL_VARIANT(path) if mp.current_process().name == 'MainProcess' else FaultFile(path, error)
                    with patch.object(pysam, 'VariantFile', opened), self.assertRaises(error):
                        if index == 0:
                            MODULES[index].parse_vcf(self.vcf, n_processes=1)
                        elif index == 1:
                            MODULES[index].main(self.vcf, str(out), str(splice), n_processes=1)
                        else:
                            with patch.object(sys, 'argv', ['stat', '--vcf_path', self.vcf, '--output', str(out), '--threads', '1']):
                                MODULES[index].main()
                    self.assertEqual(before, [(p.read_bytes(), p.stat().st_mtime_ns) for p in (out, splice)])

    def test_real_cli_missing_index_nonzero_without_publication(self):
        Path(self.vcf + '.tbi').unlink()
        for index in range(3):
            with self.subTest(producer=NAMES[index]):
                out = self.root / ('new' + str(index) + '.pkl')
                if index == 0:
                    args = ['--vcf_path', self.vcf, '--output', str(out), '--processes', '1',
                            '--rmc_tsv', 'unused', '--hcseeker_tsv', 'unused']
                elif index == 1:
                    args = [self.vcf, str(out), str(out) + '.splice', '1']
                else:
                    args = ['--vcf_path', self.vcf, '--output', str(out), '--threads', '1']
                result = subprocess.run([sys.executable, '-B', str(SCRIPTS / (NAMES[index] + '.py')), *args],
                                        capture_output=True, text=True, timeout=60)
                self.assertNotEqual(result.returncode, 0, result.stderr)
                self.assertIn('fetch requires an index', result.stderr)
                self.assertFalse(out.exists())

    def test_exon_missing_required_header_is_nonzero_before_publication(self):
        raw = self.root / 'tiny.vcf'
        raw.write_text(raw.read_text().replace('Gene|Feature|Consequence', 'Gene|Other|Consequence')
                       .replace('|EXON|', '|OtherExon|'))
        pysam.tabix_compress(str(raw), self.vcf, force=True)
        pysam.tabix_index(self.vcf, preset='vcf', force=True)
        baseline = os.environ.get('PRIVA_WORKER_BASELINE')
        if baseline:
            old_output = self.root / 'old-schema.pkl'
            old = subprocess.run([sys.executable, '-B', str(Path(baseline) / 'scripts/stat_exon_patho_afs.py'),
                                  '--vcf_path', self.vcf, '--output', str(old_output), '--threads', '1'],
                                 capture_output=True, text=True, timeout=60)
            self.assertEqual(old.returncode, 0)
            self.assertFalse(old_output.exists())
        for existing in (False, True):
            with self.subTest(existing_output=existing):
                out = self.root / ('schema-' + str(existing) + '.pkl')
                if existing:
                    out.write_bytes(b'previous complete output')
                before = (out.read_bytes(), out.stat().st_mtime_ns) if existing else None
                result = subprocess.run([sys.executable, '-B', str(SCRIPTS / 'stat_exon_patho_afs.py'),
                                         '--vcf_path', self.vcf, '--output', str(out), '--threads', '1'],
                                        capture_output=True, text=True, timeout=60)
                self.assertNotEqual(result.returncode, 0)
                self.assertIn('Required CSQ fields not found', result.stderr)
                if existing:
                    self.assertEqual(before, (out.read_bytes(), out.stat().st_mtime_ns))
                else:
                    self.assertFalse(out.exists())


if __name__ == '__main__':
    unittest.main()
