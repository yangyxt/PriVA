"""Real-tool cold GRC/UCSC equivalence; synthetic sequence, no VEP or models."""
import gzip
import importlib.util
import os
from pathlib import Path
import shutil
import subprocess
import tempfile
import unittest

ROOT = Path(__file__).resolve().parents[1]
READY = all(shutil.which(x) for x in ('bcftools', 'tabix')) and all(
    importlib.util.find_spec(x) for x in ('pysam', 'pandas'))


@unittest.skipUnless(READY, 'bcftools, tabix, pysam and pandas required')
class PreprocessContigMapping(unittest.TestCase):
    def setUp(self):
        import pysam
        self.tmp = tempfile.TemporaryDirectory(prefix='priva_contigs_', dir=os.environ.get('PRIVA_TEST_TMPDIR', os.environ.get('TMPDIR')))
        self.addCleanup(self.tmp.cleanup)
        self.root = Path(self.tmp.name)
        self.reference = self.root / 'synthetic.fa'
        self.reference.write_text('>chr1\n' + 'CGTACGTACG' + 'A' * 20 + 'CGTACGTACG' * 17 + '\n')
        pysam.faidx(str(self.reference))

    def preprocess(self, chrom, genotypes):
        import pysam
        folder = self.root / chrom
        folder.mkdir()
        source = folder / 'input.vcf'
        header = ('##fileformat=VCFv4.2\n##contig=<ID=' + chrom + ',length=200>\n'
                  '##INFO=<ID=CASE,Number=1,Type=String,Description="Synthetic case">\n')
        if genotypes:
            header += ('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n'
                       '##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Depth">\n')
        header += '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO' + ('\tFORMAT\tS' if genotypes else '') + '\n'
        rows = []
        for number, (pos, ref, alt) in enumerate(((3, 'TA', 'CG'), (8, 'A', 'G'), (20, 'AA', 'A'))):
            row = [chrom, str(pos), 'old_' + str(number), ref, alt, '40', 'PASS', 'CASE=case_' + str(number)]
            if genotypes:
                row += ['GT:AD', '0/1:10,10']
            rows.append('\t'.join(row) + '\n')
        source.write_text(header + ''.join(rows))
        compressed = str(source) + '.gz'
        pysam.tabix_compress(str(source), compressed)
        pysam.tabix_index(compressed, preset='vcf')
        target = folder / 'result.vcf.gz'
        result = subprocess.run(['bash', str(ROOT / 'scripts/annotation_vcf.sh'), 'preprocess_vcf',
                                 '-i', compressed, '-o', str(target), '-r', str(self.reference), '-t', '1'],
                                env=dict(os.environ, TMPDIR=str(self.root)), capture_output=True, text=True, timeout=60)
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertTrue(Path(str(target) + '.tbi').is_file())
        with gzip.open(target, 'rt') as stream:
            return [line.rstrip('\n').split('\t') for line in stream if not line.startswith('#')]

    def check_equivalence(self, genotypes):
        ucsc, grc = self.preprocess('chr1', genotypes), self.preprocess('1', genotypes)
        self.assertEqual(ucsc, grc)  # Includes all IDs, INFO, GT/AD and record order.
        self.assertEqual([(r[0], r[1], r[3], r[4]) for r in grc],
                         [('chr1', '3', 'T', 'C'), ('chr1', '4', 'A', 'G'),
                          ('chr1', '8', 'A', 'G'), ('chr1', '10', 'GA', 'G')])
        self.assertEqual([r[2] for r in grc], ['chr1:3:T:C', 'chr1:4:A:G', 'chr1:8:A:G', 'chr1:10:GA:G'])
        if genotypes:
            self.assertTrue(all(r[8:] == ['GT:AD', '0/1:10,10'] for r in grc))
        else:
            self.assertTrue(all(len(r) == 8 for r in grc))

    def test_genotyped_cold_grc_equals_ucsc(self):
        self.check_equivalence(True)

    def test_sites_without_format_cold_grc_equals_ucsc(self):
        self.check_equivalence(False)


if __name__ == '__main__':
    unittest.main()
