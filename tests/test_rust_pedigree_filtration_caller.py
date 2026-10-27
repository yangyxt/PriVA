"""Actual shell/task seams with explicit unrelated-tool stubs, no family data."""
import argparse
import contextlib
import importlib.util
import json
import os
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest
from unittest.mock import patch

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'scripts'))
spec = importlib.util.spec_from_file_location('pedigree_test_priva_task', ROOT / 'bin/priva_task.py')
task = importlib.util.module_from_spec(spec)
spec.loader.exec_module(task)


class PedigreeCallerControls(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory()
        self.addCleanup(self.temporary.cleanup)
        self.root = Path(self.temporary.name)
        self.binary = self.root / 'qualified-executable'
        self.binary.write_bytes(b'#!/bin/sh\nexit 0\n')
        self.binary.chmod(0o700)
        self.input = self.root / 'input.vcf.gz'
        self.input.write_text('synthetic input')
        self.output = self.root / 'output.vcf.gz'
        self.ped = self.root / 'input.ped'
        self.ped.write_text('F\tp\t0\t0\t1\t2\n')
        self.events = self.root / 'events.txt'
        self.argv = self.root / 'argv.bin'

    def shell(self, *, selected='', warm=False, python_rc=0, tabix_rc=0, display_rc=0, main=False, ranking_only=False):
        if warm:
            self.output.write_text('synthetic old output')
            os.utime(self.input, ns=(1_000_000_000, 1_000_000_000))
            os.utime(self.output, ns=(2_000_000_000, 2_000_000_000))
        script = r'''
source "$SOURCE"
log() { :; }
check_vcf_validity() { if [[ "$WARM" == 1 ]]; then printf 1; fi; }
count_vcf_records() { printf 1; }
python3() { printf 'python\n' >> "$EVENTS"; printf '%s\0' "$@" > "$ARGV"; return "$PYTHON_RC"; }
tabix() { printf 'tabix\n' >> "$EVENTS"; return "$TABIX_RC"; }
display_vcf() { printf 'display\n' >> "$EVENTS"; return "$DISPLAY_RC"; }
if [[ "$MAIN" == 1 ]]; then
 read_yaml() {
  case "$2" in
   ped_file) printf '%s\n' "$PED";; threads_per_fam) echo 1;;
   priva_pedigree_binary) printf '%s\n' "$SELECTED";;
   priva_ranking_binary) printf '%s\n' "$RANKING";;
   output_dir) printf '%s\n' "$DIR";; input_vcf) printf '%s\n' "$INPUT";;
   *) echo null;;
  esac
 }
 extract_fam_vcf() { printf 'extract\n' >> "$EVENTS"; printf fixture > "${INPUT/.vcf*/.F.vcf.gz}"; }
 filter_af() { printf 'af\n' >> "$EVENTS"; }
 main_filtration -v "$INPUT" -f F -c "$CONFIG"
else
 args=(-i "$INPUT" -p "$PED" -f F -o "$OUTPUT")
 [[ -z "$SELECTED" ]] || args+=(-b "$SELECTED")
 filter_allele_based_on_pedigree_data "${args[@]}"
fi
'''
        source = ROOT / 'scripts/filtration_vcf_per_fam.sh'
        env = dict(os.environ, SOURCE=str(source), INPUT=str(self.input), OUTPUT=str(self.output),
            PED=str(self.ped), SELECTED=str(selected), RANKING=str(self.binary) if ranking_only else '',
            WARM=str(int(warm)), PYTHON_RC=str(python_rc), TABIX_RC=str(tabix_rc), DISPLAY_RC=str(display_rc),
            MAIN=str(int(main)), CONFIG=str(self.root / 'config.yaml'), DIR=str(self.root),
            EVENTS=str(self.events), ARGV=str(self.argv))
        result = subprocess.run(['bash', '-c', script, str(source)], env=env, capture_output=True, text=True, timeout=15)
        events = self.events.read_text().splitlines() if self.events.exists() else []
        argv = self.argv.read_bytes().split(b'\0')[:-1] if self.argv.exists() else []
        return result, events, [value.decode() for value in argv]

    def test_default_actual_shell_argv_and_tool_order_are_unchanged(self):
        result, events, argv = self.shell()
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertEqual(events, ['python', 'tabix', 'display'])
        self.assertEqual(argv, [str(ROOT / 'scripts/pedigree_filtration_per_fam.py'), '-v', str(self.input),
                               '-p', str(self.ped), '-f', 'F', '-o', str(self.output)])

    def test_explicit_selection_appends_only_the_native_argument(self):
        result, events, argv = self.shell(selected=self.binary)
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertEqual(events, ['python', 'tabix', 'display'])
        self.assertEqual(argv[-2:], ['--pedigree-binary', str(self.binary)])

    def test_warm_guard_skips_both_backends_and_downstream_tools(self):
        result, events, argv = self.shell(selected=self.root / 'missing-binary', warm=True)
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertEqual((events, argv), ([], []))

    def test_cold_invalid_selection_fails_without_starting_python(self):
        result, events, _ = self.shell(selected=self.root / 'missing-binary')
        self.assertEqual(result.returncode, 1)
        self.assertEqual(events, [])

    def test_python_failure_stops_tabix_and_display(self):
        result, events, _ = self.shell(selected=self.binary, python_rc=7)
        self.assertEqual(result.returncode, 1)
        self.assertEqual(events, ['python'])

    def test_tabix_failure_stops_display(self):
        result, events, _ = self.shell(tabix_rc=9)
        self.assertEqual(result.returncode, 1)
        self.assertEqual(events, ['python', 'tabix'])

    def test_display_failure_is_not_success(self):
        result, events, _ = self.shell(display_rc=11)
        self.assertEqual(result.returncode, 1)
        self.assertEqual(events, ['python', 'tabix', 'display'])

    def test_actual_main_forwards_selected_key_and_preserves_stage_order(self):
        result, events, argv = self.shell(main=True, selected=self.binary)
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertEqual(events, ['extract', 'python', 'tabix', 'display', 'af'])
        self.assertEqual(argv[-2:], ['--pedigree-binary', str(self.binary)])

    def test_ranking_selection_does_not_activate_pedigree(self):
        result, events, argv = self.shell(main=True, ranking_only=True)
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertNotIn('--pedigree-binary', argv)
        self.assertEqual(events, ['extract', 'python', 'tabix', 'display', 'af'])

    def test_actual_main_stops_before_AF_on_pedigree_failure(self):
        result, events, _ = self.shell(main=True, selected=self.binary, python_rc=7)
        self.assertEqual(result.returncode, 1)
        self.assertEqual(events, ['extract', 'python'])

    def test_selected_resource_always_has_SHA_even_without_full_checksums(self):
        record = task.inventory({'priva_pedigree_binary': self.binary}, checksum=False)['priva_pedigree_binary']
        self.assertEqual(len(record['sha256']), 64)
        self.assertEqual(record['executable_mode'], 0o700)

    def manifest(self):
        config = self.root / 'config.yaml'
        config.write_text('assembly: hg38\n')
        output = self.root / 'manifest.json'
        args = argparse.Namespace(config=str(config), legacy_root=str(ROOT), mode='annotated',
            vcf_hub_mode='off', cadd_hub_mode='copy', checksum=False, output=str(output))
        with patch.object(task, 'runtime_resources', return_value={'priva_pedigree_binary': self.binary}):
            task.build_manifest(args)
        return config, json.loads(output.read_bytes())

    def test_manifest_catches_same_size_same_mtime_executable_change(self):
        config, manifest = self.manifest()
        self.assertEqual(manifest['identity_kind'], 'size-mtime-inventory+pedigree-sha256')
        before = self.binary.stat()
        self.binary.write_bytes(b'#!/bin/sh\nexit 1\n')
        os.utime(self.binary, ns=(before.st_atime_ns, before.st_mtime_ns))
        with patch.object(task, 'runtime_resources', return_value={'priva_pedigree_binary': self.binary}):
            with self.assertRaisesRegex(task.ContractError, 'Resource changed.*priva_pedigree_binary'):
                task.verify_manifest(manifest, config, ROOT, 'annotated', 'off', 'copy')

    def test_resource_hash_race_is_rejected(self):
        original = task.sha256
        def changed(path):
            result = original(path)
            self.binary.write_bytes(b'#!/bin/sh\nexit 1\n')
            return result
        with patch.object(task, 'sha256', side_effect=changed):
            with self.assertRaisesRegex(task.ContractError, 'changed during identity capture'):
                task.executable_record(self.binary, 'priva_pedigree_binary')

    def test_private_config_forwards_exact_binary_without_copy_or_activation(self):
        import yaml
        bundle = self.root / 'bundle'
        bundle.mkdir()
        (bundle / 'family.ped').write_bytes(self.ped.read_bytes())
        config = dict(assembly='hg38', threads=1, threads_per_fam=1,
                      priva_pedigree_binary=str(self.binary), priva_ranking_binary='')
        (bundle / 'legacy.yaml').write_text(yaml.safe_dump(config))
        before = self.binary.stat()
        previous = Path.cwd()
        try:
            os.chdir(self.root)
            path, _ = task.private_config(bundle, dict(base='small', mode='annotated', assembly='hg38'), 'filtration', 1)
            actual = yaml.safe_load(path.read_text())
        finally:
            os.chdir(previous)
        self.assertEqual(actual['priva_pedigree_binary'], str(self.binary))
        self.assertEqual(actual['priva_ranking_binary'], '')
        self.assertEqual(self.binary.stat().st_ino, before.st_ino)
        self.assertEqual(self.binary.stat().st_mtime_ns, before.st_mtime_ns)
        self.assertFalse((self.root / 'result/private_inputs/priva_pedigree_binary').exists())


if __name__ == '__main__':
    unittest.main()
