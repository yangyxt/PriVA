"""Fixture-only checks of PriVA's upstream CADD.sh call (Calculate_CADD).

A stub CADD.sh records its arguments; no CADD environment or model is run.
"""
import os
import gzip
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile
import unittest
from test_cadd_conda import write_conda

ROOT = Path(__file__).resolve().parents[1]

CADD_STUB = '''#!/usr/bin/env bash
printf '%s\\n' "$@" > "$FIXTURE/cadd-argv.txt"
printf '%s' "${PRIVA_CADD_MMSPLICE_EXON_CACHE-unset}" > "$FIXTURE/cadd-exon-tables.txt"
while getopts ':o:g:c:r:amqpd' option; do [[ $option == o ]] && out=$OPTARG; done
printf '## CADD fixture\\n#Chrom\\tPos\\tRef\\tAlt\\tRawScore\\tPHRED\\n1\\t10\\tA\\tG\\t0.5\\t12.0\\n' | gzip -c > "$out"
'''

SHELL = '''source "$ROOT/scripts/annotation_vcf.sh" > /dev/null 2>&1
SCRIPT_DIR="$FIXTURE/scripts"
read_yaml() { case "$2" in threads) printf 2;; cadd_base_dir) printf '%s/cadd' "$FIXTURE";;
    assembly) printf hg38;; tmp_dir) printf '%s/tmp' "$FIXTURE";; esac; }
check_vcf_validity() { return 0; }
display_table() { return 0; }
liftover_from_ucsc_to_GRCh() { cp -- "$1" "$3"; }
Calculate_CADD "$FIXTURE/input.vcf.gz" "$FIXTURE/config.yaml" "$FIXTURE/output.cadd.tsv"
'''


class UpstreamCaddRouteTests(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory(prefix='cadd-upstream-fixture-')
        self.root = Path(self.temporary.name)
        for name in ('cadd', 'tmp', 'scripts', 'bin'):
            (self.root / name).mkdir()
        script = self.root / 'cadd/CADD.sh'
        script.write_text(CADD_STUB)
        script.chmod(0o755)
        self.exon_tables = self.root / 'cadd/data/annotations/GRCh38_v1.7/mmsplice'
        self.exon_tables.mkdir(parents=True)
        (self.root / 'input.vcf.gz').write_bytes(b'fixture input')
        (self.root / 'config.yaml').write_text('assembly: hg38\n')
        (self.root / 'bin/python').symlink_to(sys.executable)
        (self.root / 'bin/mawk').symlink_to(shutil.which('mawk') or shutil.which('awk'))

    def tearDown(self):
        self.temporary.cleanup()

    def tool(self, name):
        path = self.root / 'bin' / name
        path.write_text('#!/usr/bin/env bash\nexit 0\n')
        path.chmod(0o755)

    def run_route(self, shell=SHELL):
        path = ':'.join([str(self.root / 'bin'), '/usr/bin', '/bin'])
        if shutil.which('conda', path='/usr/bin:/bin'):
            self.skipTest('a host conda outside the fixture would select the conda route')
        env = dict(os.environ, ROOT=str(ROOT), FIXTURE=str(self.root), PATH=path, CONDA_EXE='')
        result = subprocess.run(['bash', '-c', shell], env=env, capture_output=True, text=True, timeout=60)
        argv = (self.root / 'cadd-argv.txt').read_text().splitlines() if (self.root / 'cadd-argv.txt').exists() else []
        return result, argv

    def test_conda_host_keeps_original_conda_only_command(self):
        write_conda(self.root)
        result, argv = self.run_route()
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertEqual(argv, ['-c', '2', '-a', '-p', '-m', '-d', '-g', 'GRCh38', '-o',
                                str(self.root / 'output.cadd.tsv.gz'), str(self.root / 'input.nochr.vcf.gz')])
        # PriVA's MMSplice copy saves and reloads its exon table beside CADD's GTF.
        self.assertEqual((self.root / 'cadd-exon-tables.txt').read_text(), str(self.exon_tables))
        # The consumer receives the plain table without CADD's ## header lines.
        self.assertEqual((self.root / 'output.cadd.tsv').read_text().splitlines()[0],
                         '#Chrom\tPos\tRef\tAlt\tRawScore\tPHRED')

    def test_missing_conda_stops_before_cadd_runs(self):
        result, argv = self.run_route()
        self.assertNotEqual(result.returncode, 0)
        self.assertIn('needs a conda executable', result.stderr)
        self.assertEqual(argv, [])
        self.assertFalse((self.root / 'output.cadd.tsv').exists())

    def test_runtime_finds_frontend_without_activating_its_python(self):
        write_conda(self.root / 'conda_frontend')
        result, argv = self.run_route()
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertIn('-m', argv)

    def test_installer_upstream_check_uses_same_frontend_and_original_arguments(self):
        write_conda(self.root / 'conda_frontend')
        tests = self.root / 'cadd/test'
        tests.mkdir()
        (tests / 'input.vcf.gz').write_bytes(b'fixture input')
        (tests / 'output_inclAnno_GRCh38.tsv.gz').write_bytes(gzip.compress(
            b'## CADD fixture\n#Chrom\tPos\tRef\tAlt\tRawScore\tPHRED\n1\t10\tA\tG\t0.5\t12.0\n'))
        shell = SHELL.split('Calculate_CADD ', 1)[0].replace('/scripts/annotation_vcf.sh', '/scripts/install_utils.sh')
        result, argv = self.run_route(shell + 'CADD_upstream_check "$FIXTURE/config.yaml"\n')
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertEqual(argv[:8], ['-c', '2', '-a', '-m', '-g', 'GRCh38', '-o', argv[7]])
        self.assertTrue(argv[7].startswith(str(self.root / 'tmp/cadd_upstream_check.')))
        self.assertEqual(argv[8], str(tests / 'input.vcf.gz'))
        self.assertIn('reproduced its bundled GRCh38 test PHRED scores', result.stderr)


if __name__ == '__main__':
    unittest.main()
