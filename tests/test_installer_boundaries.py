"""Guard installer calibration and liftover before scientific tools run."""
import os
from pathlib import Path
import subprocess
import tempfile
import unittest

ROOT = Path(__file__).resolve().parents[1]


class InstallerBoundaryTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory(dir=ROOT / 'tests', prefix='installer_boundary_')
        self.addCleanup(self.temp.cleanup)
        self.folder = Path(self.temp.name)

    def test_basic_vep_refuses_mixed_spliceai_versions_before_vep(self):
        cache = self.folder / 'scores.vcf'
        cache.write_text('##fileformat=VCFv4.2\n##SpliceAI2Version=2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n')
        v1 = self.folder / 'v1_scores.vcf'
        v1.write_text('##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n')
        config = self.folder / 'config.yaml'
        config.write_text('assembly: hg38\n')
        result = subprocess.run(['bash', '-c', """
            source "$1" >/dev/null
            vep() { printf 'UNEXPECTED_VEP\n'; return 99; }
            basic_vep_annotation -i unused.vcf -f "$2" -s "$3" -d "$4"
        """, 'test', str(ROOT / 'scripts/install_utils.sh'), str(config), str(cache), str(v1)], text=True, capture_output=True)
        self.assertNotEqual(result.returncode, 0)
        self.assertIn('different versions', result.stderr)
        self.assertNotIn('UNEXPECTED_VEP', result.stdout)

    def test_failed_bgz_conversion_stops_before_crossmap(self):
        for stage in ('view', 'index'):
            with self.subTest(stage=stage):
                result = subprocess.run(['bash', '-c', """
                    load() { source "$1"; }; load "$1"
                    # Failure selector is passed through the environment.
                    bcftools() { [[ $1 != "$FAIL_STAGE" ]]; }
                    display_vcf() { return 0; }
                    CrossMap() { printf 'UNEXPECTED_CROSSMAP\n'; return 99; }
                    crossmap_liftover_hg382hg19 -i input.vcf.bgz -c chain -f ref.fa
                """, 'test', str(ROOT / 'scripts/common_bash_utils.sh')],
                    env=dict(os.environ, FAIL_STAGE=stage), text=True, capture_output=True)
                self.assertNotEqual(result.returncode, 0)
                self.assertIn('Failed to prepare BGZF gnomAD input', result.stderr)
                self.assertNotIn('UNEXPECTED_CROSSMAP', result.stdout)


if __name__ == '__main__':
    unittest.main()
