"""Two-path installation, atomic registration, and verified receipt reuse."""
import contextlib
import hashlib
import io
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile
import unittest
from unittest import mock

import yaml

ROOT = Path(__file__).resolve().parents[1]
INSTALLER = ROOT / 'scripts/install_utils.sh'
FROZEN_INSTALLER = ROOT / 'tests/fixtures/installer-054b494/scripts/install_utils.sh'


class SpliceAIInstallerInputs(unittest.TestCase):
    """Shell prompt/staging contract; shared-helper format tests live separately."""
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory(prefix='spliceai_input_', dir=ROOT / 'tests')
        self.addCleanup(self.temporary.cleanup)
        self.folder = Path(self.temporary.name)
        scripts = self.folder / 'source/scripts'
        scripts.mkdir(parents=True)
        for name in ('install_utils.sh', 'common_bash_utils.sh', 'annotation_vcf.sh'):
            shutil.copyfile(ROOT / 'scripts' / name, scripts / name)
        self.installer = scripts / 'install_utils.sh'
        self.config = self.folder / 'config.yaml'
        self.calls = self.folder / 'calls.json'
        self.prompts = self.folder / 'prompts.txt'
        self.cache = self.folder / 'cache with spaces'
        self.plugins = self.folder / 'plugins'
        self.cache.mkdir()
        self.plugins.mkdir()
        self.snv = self.folder / 'snv with spaces.vcf.gz'
        self.indel = self.folder / 'indel with spaces.vcf.gz'
        self.snv.write_bytes(b'fixture: helper owns decompressed-header validation')
        self.indel.write_bytes(b'fixture: helper owns decompressed-header validation')
        self.helper = scripts / 'spliceai_cache_input.py'
        self.helper.write_text(FAKE_INPUT_HELPER)
        self.env = dict(os.environ, TEST_CONFIG=str(self.config), TEST_FOLDER=str(self.folder),
                        TEST_CALLS=str(self.calls), TEST_PROMPTS=str(self.prompts),
                        TEST_SNV=str(self.snv), TEST_INDEL=str(self.indel))
        self.configure()

    def configure(self, assembly='hg19', **extra):
        self.config.write_text(yaml.safe_dump(dict(assembly=assembly, unrelated='keep', **extra), sort_keys=False))
        self.env['TEST_ORIGINAL_CONFIG_SHA'] = hashlib.sha256(self.config.read_bytes()).hexdigest()
        return self.config.read_bytes()

    def shell(self, body):
        return subprocess.run(['bash', '-c', 'source "$1"; ' + body, 'test', str(self.installer)],
                              env=self.env, text=True, capture_output=True)

    def invoke(self):
        return self.shell(r'''
            answers=("${TEST_CONFIRM:-yes}" "$TEST_SNV" "$TEST_INDEL")
            position=0
            read() { printf '%s\n' "$*" >> "$TEST_PROMPTS"; REPLY=${answers[$position]}; position=$((position+1)); }
            SpliceAI_install "$TEST_CONFIG" "$TEST_FOLDER/cache with spaces" "$TEST_FOLDER/plugins"
        ''')

    def test_frozen_original_oracle_hash_is_not_rewritten(self):
        self.assertEqual(hashlib.sha256(FROZEN_INSTALLER.read_bytes()).hexdigest(),
                         '80107edcc5ec845e591a87957e230fa586c2a480746c904cac3d4a84e209a81f')

    def test_fresh_assembly_guidance_preserves_three_prompts_and_stages_both_paths(self):
        for assembly in ('hg19', 'GRCh37', 'hg38', 'GRCh38'):
            with self.subTest(assembly=assembly):
                self.prompts.unlink(missing_ok=True)
                self.configure(assembly)
                result = self.invoke()
                self.assertEqual(result.returncode, 0, result.stderr)
                prompts = self.prompts.read_text().splitlines()
                self.assertEqual(len(prompts), 3)
                self.assertIn('finished downloading', prompts[0])
                self.assertIn('snv prescore file', prompts[1])
                self.assertIn('indel prescore file', prompts[2])
                calls = json.loads(self.calls.read_text())
                self.assertEqual(calls, dict(command='install', config=str(self.config),
                                            cache_dir=str(self.cache), snv=str(self.snv), indel=str(self.indel)))
                self.assertEqual(yaml.safe_load(self.config.read_bytes()), dict(assembly=assembly,
                                 unrelated='keep', spliceai_snv_prescore=str(self.snv), spliceai_indel_prescore=str(self.indel)))
                self.assertIn('SpliceAI2 v2.0' if assembly in ('hg38', 'GRCh38') else 'SpliceAI v1', result.stderr)

    def test_existing_configured_paths_are_validated_by_helper_without_prompts(self):
        self.configure('GRCh38', spliceai_snv_prescore=str(self.snv), spliceai_indel_prescore=str(self.indel))
        result = self.invoke()
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertFalse(self.prompts.exists())
        self.assertTrue(self.calls.exists())

    def test_missing_second_input_never_partially_registers_first(self):
        original = self.config.read_bytes()
        self.indel.unlink()
        result = self.invoke()
        self.assertNotEqual(result.returncode, 0)
        self.assertEqual(self.config.read_bytes(), original)
        self.assertFalse(self.calls.exists())
        self.assertEqual(len(self.prompts.read_text().splitlines()), 3)

    def test_missing_first_input_retains_original_early_failure_order(self):
        original = self.config.read_bytes()
        self.snv.unlink()
        result = self.invoke()
        self.assertNotEqual(result.returncode, 0)
        self.assertEqual(self.config.read_bytes(), original)
        self.assertFalse(self.calls.exists())
        self.assertEqual(len(self.prompts.read_text().splitlines()), 2)

    def test_declined_download_confirmation_changes_nothing(self):
        original = self.config.read_bytes()
        self.env['TEST_CONFIRM'] = 'no'
        result = self.invoke()
        self.assertNotEqual(result.returncode, 0)
        self.assertEqual(self.config.read_bytes(), original)
        self.assertFalse(self.calls.exists())
        self.assertEqual(len(self.prompts.read_text().splitlines()), 1)

    def test_helper_failure_never_registers_paths(self):
        original = self.config.read_bytes()
        self.env['TEST_HELPER_FAIL'] = '1'
        result = self.invoke()
        self.assertNotEqual(result.returncode, 0)
        self.assertEqual(self.config.read_bytes(), original)

    def test_unknown_assembly_or_missing_helper_fails_before_prompts(self):
        self.configure('hg37')
        self.assertNotEqual(self.invoke().returncode, 0)
        self.assertFalse(self.prompts.exists())
        self.configure('hg19')
        self.helper.unlink()
        self.assertNotEqual(self.invoke().returncode, 0)
        self.assertFalse(self.prompts.exists())

    def test_original_plugin_stage_order_is_retained(self):
        result = self.shell(r'''
            CONDA_PREFIX="$TEST_FOLDER/envs/priva_acmg"
            for name in UTRAnnotator_install LOEUF_install Conservation_install AlphaMissense_install AlphaMissense_anno AlphaMissense_pick_missense_constrained_segments SpliceAI_install SpliceVault_install MaveDB_install PrimateAI_install; do
                eval "$name() { printf '%s\n' '$name' >> \"\$TEST_CALLS\"; }"
            done
            VEP_plugins_install plugins caches hg38 "$TEST_CONFIG" priva_acmg
        ''')
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertEqual(self.calls.read_text().splitlines(), [
            'UTRAnnotator_install', 'LOEUF_install', 'Conservation_install',
            'AlphaMissense_install', 'AlphaMissense_anno',
            'AlphaMissense_pick_missense_constrained_segments', 'SpliceAI_install',
            'SpliceVault_install', 'MaveDB_install', 'PrimateAI_install'])
        self.assertNotIn('SpliceAI_profile_install', INSTALLER.read_text())


FAKE_INPUT_HELPER = r'''
import argparse, hashlib, json, os
from pathlib import Path
import yaml
p=argparse.ArgumentParser()
p.add_argument('command', choices=['install'])
for key in ('config','cache-dir','snv','indel'):p.add_argument('--'+key,required=True)
a=p.parse_args()
config=Path(a.config)
assert hashlib.sha256(config.read_bytes()).hexdigest()==os.environ['TEST_ORIGINAL_CONFIG_SHA'], 'Installer updated config before validating both inputs'
Path(os.environ['TEST_CALLS']).write_text(json.dumps(vars(a)))
if os.environ.get('TEST_HELPER_FAIL'):raise SystemExit(17)
settings=yaml.safe_load(config.read_bytes())
settings.update(spliceai_snv_prescore=a.snv,spliceai_indel_prescore=a.indel)
tmp=config.with_suffix('.tmp');tmp.write_text(yaml.safe_dump(settings,sort_keys=False));tmp.replace(config)
'''


FAKE_DRIVER = r'''
import argparse
import hashlib
import json
import os
from pathlib import Path
import subprocess

parser = argparse.ArgumentParser()
parser.add_argument('command', choices=['prepare'])
for flag in ('manifest', 'assembly', 'output-dir', 'receipt'):
    parser.add_argument('--' + flag, required=True)
parser.add_argument('--converter-bin')
args = parser.parse_args()
Path(os.environ['TEST_CALLS']).write_text(json.dumps(vars(args)))
mode = os.environ.get('TEST_ADAPTER_MODE', 'valid')
if mode == 'failure':
    raise SystemExit(47)
output = Path(args.output_dir)
output.mkdir(parents=True, exist_ok=True)

def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()

files = {}
for role in ('snv', 'indel'):
    vcf = output / (role + '.vcf.gz')
    raw = ('##fileformat=VCFv4.2\n##contig=<ID=chr1,length=248956422>\n'
           '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\nchr1\t10\t.\tA\tC\t.\t.\t.\n')
    with vcf.open('wb') as stream:
        subprocess.run(['bgzip', '-c'], input=raw.encode(), stdout=stream, check=True)
    subprocess.run(['tabix', '-p', 'vcf', str(vcf)], check=True)
    index = Path(str(vcf) + '.tbi')
    files[role] = dict(path=str(vcf), index_path=str(index), sha256=digest(vcf),
                       index_sha256=digest(index), bytes=vcf.stat().st_size)
receipt = dict(schema_version=1, resource='SpliceAI2', version='2.0', assembly='GRCh38',
               status='prepared_annotation_only', scientific_qualified=False,
               coverage='complete_published',
               source_revision='cd1b8b82bcef9c97ac3b3b9a24f0de3a5bbf77f9',
               source_manifest_sha256=digest(Path(args.manifest)), outputs=files)
changes = {
    'subset': ('status', 'subset_validation_only'),
    'wrong_assembly': ('assembly', 'GRCh37'),
    'wrong_version': ('version', '1.3'),
    'wrong_revision': ('source_revision', '0' * 40),
    'wrong_manifest': ('source_manifest_sha256', '0' * 64),
    'scientifically_qualified': ('scientific_qualified', True),
    'false_as_zero': ('scientific_qualified', 0),
    'schema_as_bool': ('schema_version', True),
    'partial_coverage': ('coverage', 'subset'),
}
if mode in changes:
    key, value = changes[mode]
    receipt[key] = value
elif mode == 'missing_role':
    del files['indel']
elif mode == 'extra_role':
    files['extra'] = files['snv']
elif mode == 'wrong_sha':
    files['indel']['sha256'] = '0' * 64
elif mode == 'wrong_index_sha':
    files['snv']['index_sha256'] = '0' * 64
elif mode == 'wrong_size':
    files['snv']['bytes'] += 1
elif mode == 'missing_output':
    Path(files['snv']['path']).unlink()
elif mode == 'empty_index':
    Path(files['snv']['index_path']).write_bytes(b'')
elif mode == 'wrong_index_path':
    files['snv']['index_path'] = files['indel']['index_path']
    files['snv']['index_sha256'] = files['indel']['index_sha256']
elif mode == 'duplicate_outputs':
    files['indel'] = files['snv']
elif mode == 'partial_output':
    source = Path(files['snv']['path'])
    partial = output / ('snv.partial.vcf.gz')
    source.rename(partial)
    files['snv']['path'] = str(partial)
elif mode == 'outside_output':
    source = Path(files['snv']['path'])
    outside = output.parent / 'outside.vcf.gz'
    source.rename(outside)
    files['snv']['path'] = str(outside)
elif mode == 'dotdot_output':
    source = Path(files['snv']['path'])
    outside = output.parent / 'outside.vcf.gz'
    source.rename(outside)
    files['snv']['path'] = str(output / '..' / 'outside.vcf.gz')
elif mode == 'symlink_output':
    source = Path(files['snv']['path'])
    target = output / 'target.vcf.gz'
    source.rename(target)
    source.symlink_to(target)
elif mode == 'config_changed':
    with Path(os.environ['TEST_CONFIG']).open('a') as stream:
        stream.write('concurrent_setting: preserved\n')
elif mode == 'manifest_changed':
    with Path(args.manifest).open('a') as stream:
        stream.write('\n')
if mode == 'no_receipt':
    raise SystemExit(0)
if mode == 'malformed_receipt':
    Path(args.receipt).write_text('{not JSON')
elif mode == 'nonobject_receipt':
    Path(args.receipt).write_text('[]')
else:
    Path(args.receipt).write_text(json.dumps(receipt))
'''


class SpliceAI2Preparation(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory(prefix='spliceai2_prepare_', dir=ROOT / 'tests')
        self.addCleanup(self.temporary.cleanup)
        self.folder = Path(self.temporary.name)
        scripts = self.folder / 'source/scripts'
        scripts.mkdir(parents=True)
        for name in ('install_utils.sh', 'common_bash_utils.sh', 'annotation_vcf.sh'):
            shutil.copyfile(ROOT / 'scripts' / name, scripts / name)
        self.installer = scripts / 'install_utils.sh'
        self.driver = scripts / 'spliceai2_to_vep_vcf.py'
        self.driver.write_text(FAKE_DRIVER)
        self.config = self.folder / 'config.yaml'
        self.config.write_text('# Keep original settings and comments.\nassembly: hg38\n'
                               'spliceai_snv_prescore: /existing/v1-snv.vcf.gz\n'
                               'spliceai_indel_prescore: /existing/v1-indel.vcf.gz\n'
                               'unrelated: "ampersand & literal"\n')
        self.original = self.config.read_bytes()
        self.manifest = self.folder / 'complete.json'
        self.source = dict(status='downloads_complete', confirmed_by='synthetic fixture', assembly='GRCh38',
                           version='2.0', source_revision='cd1b8b82bcef9c97ac3b3b9a24f0de3a5bbf77f9',
                           files=[dict(kind='snv', chrom='chr1', path='/complete/chr1.tsv.gz')])
        # Driver tests separately qualify the complete 47-source manifest. This
        # route fixture mocks that driver and validates its verified receipt.
        self.manifest.write_text(json.dumps(self.source))
        self.output = self.folder / 'output with spaces'
        self.calls = self.folder / 'calls.json'
        self.env = dict(os.environ, TEST_CONFIG=str(self.config), TEST_CALLS=str(self.calls))

    def invoke(self, mode='valid', converter='/pinned bin/spliceai2-converter'):
        self.env['TEST_ADAPTER_MODE'] = mode
        return subprocess.run(['bash', '-c', 'source "$1"; SpliceAI2_install "$2" "$3" "$4" "$5"',
                               'test', str(self.installer), str(self.config), str(self.manifest),
                               str(self.output), converter], env=self.env, text=True, capture_output=True)

    def test_verified_preparation_registers_only_original_two_keys(self):
        self.config.chmod(0o640)
        result = self.invoke()
        self.assertEqual(result.returncode, 0, result.stderr)
        after = yaml.safe_load(self.config.read_text())
        before = yaml.safe_load(self.original)
        self.assertEqual(set(after), set(before))
        self.assertEqual({key: value for key, value in after.items()
                          if key not in ('spliceai_snv_prescore', 'spliceai_indel_prescore')},
                         {key: value for key, value in before.items()
                          if key not in ('spliceai_snv_prescore', 'spliceai_indel_prescore')})
        self.assertEqual(after['spliceai_snv_prescore'], str(self.output / 'snv.vcf.gz'))
        self.assertEqual(after['spliceai_indel_prescore'], str(self.output / 'indel.vcf.gz'))
        self.assertTrue(self.config.read_bytes().startswith(b'# Keep original settings and comments.\n'))
        self.assertEqual(self.config.stat().st_mode & 0o777, 0o640)
        calls = json.loads(self.calls.read_text())
        self.assertEqual(calls, dict(command='prepare', manifest=str(self.manifest), assembly='GRCh38',
                                     output_dir=str(self.output), receipt=str(self.output / 'adapter.ready.json'),
                                     converter_bin='/pinned bin/spliceai2-converter'))
        self.assertEqual(result.stdout.count('verify SHA256'), 5)
        self.assertIn('activation remains pending', result.stdout)

    def test_default_route_omits_production_converter_argument(self):
        result = self.invoke(converter='')
        self.assertEqual(result.returncode, 0, result.stderr)
        calls = json.loads(self.calls.read_text())
        self.assertEqual(calls['command'], 'prepare')
        self.assertIsNone(calls['converter_bin'])

    def test_invalid_receipts_never_register_cache_paths(self):
        for mode in ('no_receipt', 'malformed_receipt', 'nonobject_receipt', 'subset', 'wrong_assembly', 'wrong_version', 'wrong_revision', 'wrong_manifest',
                     'scientifically_qualified', 'false_as_zero', 'schema_as_bool', 'partial_coverage',
                     'missing_role', 'extra_role', 'wrong_sha', 'wrong_index_sha', 'wrong_size',
                     'missing_output', 'empty_index', 'wrong_index_path', 'duplicate_outputs',
                     'partial_output', 'outside_output', 'dotdot_output', 'symlink_output'):
            with self.subTest(mode=mode):
                shutil.rmtree(self.output, ignore_errors=True)
                result = self.invoke(mode)
                self.assertNotEqual(result.returncode, 0, result.stdout + result.stderr)
                self.assertEqual(self.config.read_bytes(), self.original)
                self.assertNotIn('preparation recorded', result.stdout)

    def test_incomplete_or_wrong_assembly_manifest_is_rejected_before_driver(self):
        for key, value in [('status', 'downloading'), ('confirmed_by', ''), ('assembly', 'GRCh37'),
                           ('version', '1.3'), ('source_revision', '0' * 40)]:
            with self.subTest(key=key):
                source = dict(self.source, **{key: value})
                self.manifest.write_text(json.dumps(source))
                result = self.invoke()
                self.assertNotEqual(result.returncode, 0)
                self.assertFalse(self.calls.exists())
                self.assertFalse(self.output.exists())
                self.assertEqual(self.config.read_bytes(), self.original)

    def test_partial_source_path_is_rejected_before_driver(self):
        self.source['files'][0]['path'] = '/downloads/chr1.tsv.gz.partial'
        self.manifest.write_text(json.dumps(self.source))
        result = self.invoke()
        self.assertNotEqual(result.returncode, 0)
        self.assertIn('.partial', result.stderr)
        self.assertFalse(self.calls.exists())

    def test_grch37_rejection_precedes_manifest_access(self):
        self.config.write_text('assembly: hg19\n')
        self.manifest.unlink()
        result = self.invoke()
        self.assertNotEqual(result.returncode, 0)
        self.assertIn('requires hg38/GRCh38', result.stderr)
        self.assertNotIn('FileNotFoundError', result.stderr)
        self.assertFalse(self.calls.exists())

    def test_missing_driver_and_failed_driver_do_not_update_config(self):
        result = self.invoke('failure')
        self.assertNotEqual(result.returncode, 0)
        self.assertEqual(self.config.read_bytes(), self.original)
        self.driver.unlink()
        result = self.invoke()
        self.assertNotEqual(result.returncode, 0)
        self.assertIn('deployment driver is absent', result.stderr)
        self.assertEqual(self.config.read_bytes(), self.original)

    def test_changed_config_is_not_overwritten(self):
        result = self.invoke('config_changed')
        self.assertNotEqual(result.returncode, 0)
        self.assertIn('config changed during', result.stderr)
        self.assertEqual(self.config.read_bytes(), self.original + b'concurrent_setting: preserved\n')

    def test_changed_manifest_does_not_register_outputs(self):
        result = self.invoke('manifest_changed')
        self.assertNotEqual(result.returncode, 0)
        self.assertIn('source manifest changed during', result.stderr)
        self.assertEqual(self.config.read_bytes(), self.original)

    def test_output_alias_or_partial_destination_is_rejected(self):
        real = self.folder / 'real-output'
        real.mkdir()
        self.output.symlink_to(real, target_is_directory=True)
        result = self.invoke()
        self.assertNotEqual(result.returncode, 0)
        self.assertEqual(self.config.read_bytes(), self.original)
        self.output.unlink()
        self.output = self.folder / 'output.partial'
        self.calls.unlink(missing_ok=True)
        result = self.invoke()
        self.assertNotEqual(result.returncode, 0)
        self.assertFalse(self.calls.exists())

    def test_output_changed_after_checksum_is_rejected_before_config_write(self):
        result = self.invoke()
        self.assertEqual(result.returncode, 0, result.stderr)
        self.config.write_bytes(self.original)
        receipt = self.output / 'adapter.ready.json'
        data = json.loads(receipt.read_text())
        target = Path(data['outputs']['snv']['path'])
        expected = data['outputs']['snv']['sha256']
        real_sha256 = hashlib.sha256
        changed = []

        class MutatingDigest:
            def __init__(self, *args, **kwargs):
                self.digest = real_sha256(*args, **kwargs)

            def update(self, value):
                self.digest.update(value)

            def hexdigest(self):
                value = self.digest.hexdigest()
                if value == expected and not changed:
                    with target.open('ab') as stream:
                        stream.write(b'changed after checksum')
                    changed.append(True)
                return value

        marker = '    python - "${config_file}" "${manifest}" "${output_dir}" "${receipt}" "${config_sha}" "${manifest_sha}" <<\'PY\'\n'
        body = self.installer.read_text().split(marker, 1)[1].split('\nPY\n', 1)[0]
        argv = ['receipt-check', str(self.config), str(self.manifest), str(self.output), str(receipt),
                real_sha256(self.original).hexdigest(), real_sha256(self.manifest.read_bytes()).hexdigest()]
        with mock.patch.object(sys, 'argv', argv), mock.patch('hashlib.sha256', MutatingDigest):
            with contextlib.redirect_stdout(io.StringIO()):
                with self.assertRaisesRegex(SystemExit, 'verification input changed'):
                    exec(compile(body, '<SpliceAI2 receipt verifier>', 'exec'), {})
        self.assertEqual(changed, [True])
        self.assertEqual(self.config.read_bytes(), self.original)

    def test_completed_receipt_reuse_never_invokes_converter(self):
        self.assertEqual(self.invoke().returncode, 0)
        self.calls.unlink()
        self.driver.unlink()
        self.config.write_bytes(self.original)
        result = self.invoke()
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertFalse(self.calls.exists())
        self.assertIn('Reusing the completed', result.stderr)
        self.assertEqual(yaml.safe_load(self.config.read_bytes())['spliceai_snv_prescore'],
                         str(self.output / 'snv.vcf.gz'))

    def test_invalid_existing_receipt_is_rejected_without_duplicate_conversion(self):
        self.assertEqual(self.invoke().returncode, 0)
        self.calls.unlink()
        self.config.write_bytes(self.original)
        (self.output / 'adapter.ready.json').write_text('{}')
        result = self.invoke()
        self.assertNotEqual(result.returncode, 0)
        self.assertFalse(self.calls.exists())
        self.assertEqual(self.config.read_bytes(), self.original)

    def test_incomplete_existing_output_is_not_overwritten(self):
        self.output.mkdir()
        pending = self.output / 'unfinished.vcf.gz'
        pending.write_bytes(b'not complete')
        result = self.invoke()
        self.assertNotEqual(result.returncode, 0)
        self.assertFalse(self.calls.exists())
        self.assertEqual(pending.read_bytes(), b'not complete')
        self.assertEqual(self.config.read_bytes(), self.original)

    def test_building_directory_prevents_reuse_even_with_receipt(self):
        self.assertEqual(self.invoke().returncode, 0)
        self.calls.unlink()
        self.config.write_bytes(self.original)
        (self.output / '.building-active').mkdir()
        result = self.invoke()
        self.assertNotEqual(result.returncode, 0)
        self.assertFalse(self.calls.exists())
        self.assertEqual(self.config.read_bytes(), self.original)

    def test_active_install_lock_prevents_receipt_reuse(self):
        import fcntl
        self.assertEqual(self.invoke().returncode, 0)
        self.calls.unlink()
        self.config.write_bytes(self.original)
        with (self.output / '.install.lock').open('a') as lock:
            fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
            result = self.invoke()
        self.assertNotEqual(result.returncode, 0)
        self.assertIn('preparation is active', result.stderr)
        self.assertFalse(self.calls.exists())
        self.assertEqual(self.config.read_bytes(), self.original)

    def test_duplicate_original_path_key_rejected_without_rewriting_config(self):
        self.config.write_bytes(self.original + b'spliceai_snv_prescore: duplicate\n')
        original = self.config.read_bytes()
        result = self.invoke()
        self.assertNotEqual(result.returncode, 0)
        self.assertIn('Duplicate SpliceAI2 config key', result.stderr)
        self.assertEqual(self.config.read_bytes(), original)


if __name__ == '__main__':
    unittest.main()
