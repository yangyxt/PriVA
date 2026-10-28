"""Fixture-only checks of the CADD route used when no qualified profile is selected.

A stub CADD.sh records its arguments; no CADD environment or model is run.
"""
import os
from pathlib import Path
import shutil
import subprocess
import tempfile
import unittest

ROOT = Path(__file__).resolve().parents[1]
BIN = Path('/stornext/snfs170/next-gen/scratch/yangyxt/mamba/envs/priva_acmg/bin')

CADD_STUB = '''#!/usr/bin/env bash
printf '%s\\n' "$@" > "$FIXTURE/cadd-argv.txt"
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
        # The qualified-profile helper must never be reached on this route.
        (self.root / 'scripts/cadd_shared_workflow.py').write_text(
            'import pathlib, os\npathlib.Path(os.environ["FIXTURE"], "profile-helper-called").touch()\n')
        (self.root / 'input.vcf.gz').write_bytes(b'fixture input')
        (self.root / 'config.yaml').write_text('assembly: hg38\n')

    def tearDown(self):
        self.temporary.cleanup()

    def tool(self, name):
        path = self.root / 'bin' / name
        path.write_text('#!/usr/bin/env bash\nexit 0\n')
        path.chmod(0o755)

    def run_route(self):
        path = ':'.join([str(self.root / 'bin'), str(BIN), '/usr/bin', '/bin'])
        if shutil.which('conda', path=':'.join([str(BIN), '/usr/bin', '/bin'])):
            self.skipTest('a host conda outside the fixture would select the conda route')
        env = dict(os.environ, ROOT=str(ROOT), FIXTURE=str(self.root), PATH=path)
        result = subprocess.run(['bash', '-c', SHELL], env=env, capture_output=True, text=True, timeout=60)
        argv = (self.root / 'cadd-argv.txt').read_text().splitlines() if (self.root / 'cadd-argv.txt').exists() else []
        return result, argv

    def test_conda_host_keeps_original_conda_only_command(self):
        self.tool('conda')
        result, argv = self.run_route()
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertEqual(argv, ['-c', '2', '-a', '-p', '-m', '-d', '-g', 'GRCh38', '-o',
                                str(self.root / 'output.cadd.tsv.gz'), str(self.root / 'input.nochr.vcf.gz')])
        self.assertFalse((self.root / 'profile-helper-called').exists())
        # The consumer receives the plain table without CADD's ## header lines.
        self.assertEqual((self.root / 'output.cadd.tsv').read_text().splitlines()[0],
                         '#Chrom\tPos\tRef\tAlt\tRawScore\tPHRED')

    def test_missing_conda_stops_before_cadd_runs(self):
        result, argv = self.run_route()
        self.assertNotEqual(result.returncode, 0)
        self.assertIn('needs a conda executable', result.stderr)
        self.assertEqual(argv, [])
        self.assertFalse((self.root / 'output.cadd.tsv').exists())
        self.assertFalse((self.root / 'profile-helper-called').exists())


if __name__ == '__main__':
    unittest.main()
