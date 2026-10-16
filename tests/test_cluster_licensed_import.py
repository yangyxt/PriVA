"""Preflight barriers must run before reading any licensed cache file."""
import contextlib
import hashlib
import io
import json
import os
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest
from unittest import mock

import yaml

ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / 'scripts/import_cluster_licensed.sbatch'
BODY = SCRIPT.read_text().split("<<'PY'\n", 1)[1].split('\nPY\n', 1)[0]

# Exact bytes audited from this historical installer revision. CI uses a shallow
# checkout and must neither require the worker Git object nor silently change it.
REFERENCE_REVISION = '054b4944055054148005016192647a8ed0384eef'
REFERENCE_SHA256 = {
    'install_utils.sh': '80107edcc5ec845e591a87957e230fa586c2a480746c904cac3d4a84e209a81f',
    'common_bash_utils.sh': '4cd2cffed0d6a46d1280eac97b918c1dd6c717c881dc41b0a819fbecaed439b9',
    'annotation_vcf.sh': '6b7ad6fd8eb118841cf36819193f7ce2ee14363432c66f0da75d7788c2d8d280',
    'deploy_cluster_resources.py': 'c5e345f6f3af31aff1e2a724e248f0f32a1085afe7b347183022196b359ea6cf',
}


def copy_frozen_scripts(destination, checkout=ROOT):
    contents = {}
    for name, expected in REFERENCE_SHA256.items():
        content = (checkout / 'scripts' / name).read_bytes()
        actual = hashlib.sha256(content).hexdigest()
        if actual != expected:
            raise AssertionError(
                f'Frozen installer reference mismatch: {name}: expected {expected}, got {actual}. '
                f'Reference is {REFERENCE_REVISION}; review and preserve its source bytes '
                'before changing this oracle. No current-source fallback is permitted.')
        contents[name] = content
    destination.mkdir(parents=True)
    for name, content in contents.items():
        (destination / name).write_bytes(content)


class LicensedImportPreflight(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory(prefix='licensed_import_', dir=ROOT / 'tests')
        self.addCleanup(self.temporary.cleanup)
        self.folder = Path(self.temporary.name)
        self.template = self.folder / 'template.yaml'
        self.template.write_text('assembly: hg19\n')
        self.config = self.folder / 'separate.yaml'
        self.confirmation = self.folder / 'confirmed.json'

    def run_import(self, confirmation, resource_set=None):
        self.confirmation.write_text(json.dumps(confirmation))
        return subprocess.run(
            [sys.executable, '-c', BODY, str(ROOT), '054b4944055054148005016192647a8ed0384eef',
             str(self.template), str(self.config), str(self.confirmation), str(self.folder)] +
            ([resource_set] if resource_set is not None else []),
            text=True, capture_output=True)

    def test_incomplete_acquisition_stops_before_template_or_cache_access(self):
        self.template.unlink()
        result = self.run_import({'status': 'downloading', 'confirmed_by': 'owner'})
        self.assertNotEqual(result.returncode, 0)
        self.assertIn('Explicit acquisition-owner completion confirmation', result.stderr)
        self.assertNotIn('FileNotFoundError', result.stderr)
        self.assertFalse(self.config.exists())

    def test_missing_owner_confirmation_is_rejected(self):
        result = self.run_import({'status': 'downloads_complete'})
        self.assertNotEqual(result.returncode, 0)
        self.assertIn('Explicit acquisition-owner completion confirmation', result.stderr)
        self.assertFalse(self.config.exists())

    def test_wrong_assembly_stops_before_cache_access(self):
        result = self.run_import({'status': 'downloads_complete', 'confirmed_by': 'owner',
                                  'assembly': 'hg38', 'files': {}})
        self.assertNotEqual(result.returncode, 0)
        self.assertIn('Confirmed files belong to another assembly', result.stderr)
        self.assertFalse(self.config.exists())

    def test_missing_resource_roles_stop_before_cache_access(self):
        result = self.run_import({'status': 'downloads_complete', 'confirmed_by': 'owner',
                                  'assembly': 'GRCh37', 'files': {}})
        self.assertNotEqual(result.returncode, 0)
        self.assertIn('All three matching-assembly resources are required', result.stderr)
        self.assertFalse(self.config.exists())

    def test_existing_config_is_preserved(self):
        self.config.write_text('existing configuration\n')
        result = self.run_import({'status': 'downloads_complete', 'confirmed_by': 'owner',
                                  'assembly': 'hg19', 'files': {}})
        self.assertNotEqual(result.returncode, 0)
        self.assertIn('Refusing to replace an existing import config', result.stderr)
        self.assertEqual(self.config.read_text(), 'existing configuration\n')

    def test_spliceai_selection_requires_both_roles(self):
        result = self.run_import({'status': 'downloads_complete', 'confirmed_by': 'owner',
                                  'assembly': 'hg19', 'files': {'spliceai_snv': {}}}, 'spliceai')
        self.assertNotEqual(result.returncode, 0)
        self.assertIn('Exactly the selected matching-assembly resource roles', result.stderr)
        self.assertNotIn('KeyError', result.stderr)
        self.assertFalse(self.config.exists())

    def test_spliceai_subset_is_not_the_implicit_default(self):
        result = self.run_import({'status': 'downloads_complete', 'confirmed_by': 'owner',
                                  'assembly': 'hg19',
                                  'files': {'spliceai_snv': {}, 'spliceai_indel': {}}})
        self.assertNotEqual(result.returncode, 0)
        self.assertIn('All three matching-assembly resources are required', result.stderr)
        self.assertFalse(self.config.exists())

    def test_unselected_confirmation_role_rejected_before_cache_access(self):
        result = self.run_import({'status': 'downloads_complete', 'confirmed_by': 'owner',
                                  'assembly': 'hg19', 'files': {
                                      'spliceai_snv': {}, 'spliceai_indel': {}, 'primateai': {}}}, 'spliceai')
        self.assertNotEqual(result.returncode, 0)
        self.assertIn('Exactly the selected matching-assembly resource roles', result.stderr)
        self.assertFalse(self.config.exists())

    def test_explicit_confirmation_selection_must_match(self):
        result = self.run_import({'status': 'downloads_complete', 'confirmed_by': 'owner',
                                  'resources': ['primateai']}, 'spliceai')
        self.assertNotEqual(result.returncode, 0)
        self.assertIn('Confirmed resource selection does not match', result.stderr)
        self.assertFalse(self.config.exists())

    def test_unknown_selection_is_rejected(self):
        result = self.run_import({}, 'spliceai2')
        self.assertNotEqual(result.returncode, 0)
        self.assertIn('Unknown resource set', result.stderr)
        self.assertFalse(self.config.exists())


class LicensedImportFixtures(unittest.TestCase):
    """Tiny BGZF fixtures exercise actual frozen installer prompts; no licensed files."""
    REVISION = REFERENCE_REVISION

    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory(prefix='licensed_fixture_', dir=ROOT / 'tests')
        self.addCleanup(self.temporary.cleanup)
        self.folder = Path(self.temporary.name)
        self.source = self.folder / 'frozen'
        scripts = self.source / 'scripts'
        copy_frozen_scripts(scripts)
        self.cache = self.folder / 'caches'
        spliceai = self.cache / 'SpliceAI_v1'
        spliceai.mkdir(parents=True)
        plugin_cache = self.folder / 'plugin_cache'
        plugins = self.folder / 'plugins'
        plugin_cache.mkdir()
        plugins.mkdir()
        self.template = self.folder / 'template.yaml'
        self.settings = {'assembly': 'hg19', 'vep_plugins_cachedir': str(plugin_cache),
                         'vep_plugins_dir': str(plugins), 'spliceai_snv_prescore': 'old_snv',
                         'spliceai_indel_prescore': 'old_indel', 'primateai_prescore': '/unverified/do-not-read.tsv.gz',
                         'unrelated_setting': {'keep': ['original', 17]}}
        self.template.write_text(yaml.safe_dump(self.settings))
        self.config = self.folder / 'separate.yaml'
        self.confirmation = self.folder / 'confirmed.json'
        self.confirmed = {'status': 'downloads_complete', 'confirmed_by': 'fixture_owner',
                          'assembly': 'GRCh37', 'resources': ['spliceai'], 'files': {}}
        for variant in ('snv', 'indel'):
            path = spliceai / f'spliceai_scores.raw.{variant}.hg19.vcf.gz'
            self.bgzip(path, '##fileformat=VCFv4.2\n##contig=<ID=1>\n'
                       '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n1\t1\t.\tA\tG\t.\tPASS\t.\n')
            subprocess.run(['tabix', '-p', 'vcf', str(path)], check=True)
            self.confirmed['files']['spliceai_' + variant] = {'path': str(path), 'bytes': path.stat().st_size,
                                                             'sha256': self.sha(path)}
        namespace = {'__name__': 'licensed_import_fixture'}
        exec(compile(BODY, str(SCRIPT), 'exec'), namespace)
        self.import_resources = namespace['import_resources']
        self.stdout = io.StringIO()
        self.addCleanup(mock.patch.stopall)
        mock.patch.dict(os.environ, {'SLURM_JOB_ID': 'fixture'}).start()

    @staticmethod
    def sha(path):
        return hashlib.sha256(path.read_bytes()).hexdigest()

    @staticmethod
    def bgzip(path, content):
        with path.open('wb') as stream:
            subprocess.run(['bgzip', '-c'], input=content.encode(), stdout=stream, check=True)

    def call_import(self, selection='spliceai'):
        self.confirmation.write_text(json.dumps(self.confirmed))
        with contextlib.redirect_stdout(self.stdout):
            return self.import_resources(self.source, self.REVISION, self.template, self.config,
                                         self.confirmation, self.folder, selection, cache_root=self.cache)

    def assert_rejected(self, message, selection='spliceai'):
        with self.assertRaisesRegex(RuntimeError, message):
            self.call_import(selection)
        self.assertFalse(self.config.exists())
        self.assertEqual(list(self.folder.glob('*.ready.json')), [])

    def test_spliceai_only_uses_frozen_installer_and_preserves_other_keys_and_bytes(self):
        files = list((self.cache / 'SpliceAI_v1').iterdir())
        before = {p: (self.sha(p), p.stat().st_mtime_ns) for p in files}
        template_bytes = self.template.read_bytes()
        receipt_path = self.call_import()
        self.assertEqual(receipt_path.name, 'licensed-caches-hg19-spliceai.ready.json')
        receipt = json.loads(receipt_path.read_text())
        self.assertEqual(receipt['status'], 'selected_licensed_inputs_imported')
        self.assertEqual(receipt['selected_resources'], ['spliceai'])
        self.assertEqual(receipt['resource_status'], {'spliceai': 'imported', 'primateai': 'not_selected_unverified'})
        self.assertNotIn('formatted_primateai', receipt)
        self.assertEqual([cmd[4] for cmd in receipt['commands']], ['SpliceAI_install'])
        self.assertEqual(receipt['commands'][0][3], str(self.source / 'scripts/install_utils.sh'))
        expected = dict(self.settings)
        for role, item in self.confirmed['files'].items():
            expected[role + '_prescore'] = item['path']
        self.assertEqual(yaml.safe_load(self.config.read_text()), expected)
        self.assertEqual(self.template.read_bytes(), template_bytes)
        self.assertEqual({p: (self.sha(p), p.stat().st_mtime_ns) for p in files}, before)
        self.assertFalse((self.cache / 'PrimateAI_v0_2').exists())
        self.assertIn('SHA256 start', self.stdout.getvalue())
        self.assertIn('gzip integrity check start', self.stdout.getvalue())
        self.assertIn('indexed chromosome/query checks start', self.stdout.getvalue())

    def add_primate(self):
        root = self.cache / 'PrimateAI_v0_2'
        root.mkdir()
        path = root / 'PrimateAI_scores_v0.2.sorted.tsv.bgz'
        self.bgzip(path, '#chrom\tpos\tscore\n1\t1\t0.5\n')
        subprocess.run(['tabix', '-s', '1', '-b', '2', '-e', '2', str(path)], check=True)
        self.confirmed['files']['primateai'] = {'path': str(path)}

    def test_combined_default_keeps_original_receipt_and_function_order(self):
        self.add_primate()
        del self.confirmed['resources']  # Original confirmation schema remains valid.
        self.confirmation.write_text(json.dumps(self.confirmed))
        with contextlib.redirect_stdout(self.stdout):
            receipt_path = self.import_resources(self.source, self.REVISION, self.template, self.config,
                                                 self.confirmation, self.folder, cache_root=self.cache)
        receipt = json.loads(receipt_path.read_text())
        self.assertEqual(receipt_path.name, 'licensed-caches-hg19.ready.json')
        self.assertEqual(receipt['status'], 'licensed_inputs_imported')
        self.assertEqual([cmd[4] for cmd in receipt['commands']], ['SpliceAI_install', 'PrimateAI_install'])
        self.assertEqual(receipt['resource_status'], {'spliceai': 'imported', 'primateai': 'imported'})
        self.assertEqual(len(receipt['registered_paths']), 3)
        self.assertIn('formatted_primateai_sha256', receipt)

    def test_primate_only_leaves_spliceai_keys_unverified_and_unchanged(self):
        self.add_primate()
        self.confirmed['resources'] = ['primateai']
        self.confirmed['files'] = {'primateai': self.confirmed['files']['primateai']}
        receipt = json.loads(self.call_import('primateai').read_text())
        self.assertEqual([cmd[4] for cmd in receipt['commands']], ['PrimateAI_install'])
        self.assertEqual(receipt['resource_status']['spliceai'], 'not_selected_unverified')
        result = yaml.safe_load(self.config.read_text())
        self.assertEqual(result['spliceai_snv_prescore'], 'old_snv')
        self.assertEqual(result['spliceai_indel_prescore'], 'old_indel')

    def test_reference_mismatch_stops_before_creating_a_fixture_snapshot(self):
        script = self.source / 'scripts/install_utils.sh'
        script.write_bytes(script.read_bytes() + b'\n# unexpected revision\n')
        destination = self.folder / 'rejected_snapshot/scripts'
        with self.assertRaisesRegex(AssertionError, 'Frozen installer reference mismatch'):
            copy_frozen_scripts(destination, checkout=self.source)
        self.assertFalse(destination.exists())

    def test_resource_path_must_be_absolute(self):
        self.confirmed['files']['spliceai_indel']['path'] = 'relative.vcf.gz'
        self.assert_rejected('paths must be absolute')

    def test_resource_path_outside_expected_root_is_rejected(self):
        self.confirmed['files']['spliceai_indel']['path'] = str(self.folder / 'elsewhere.vcf.gz')
        self.assert_rejected('Unexpected shared-cache root')

    def test_resource_role_cannot_point_at_the_other_variant_type(self):
        self.confirmed['files']['spliceai_indel']['path'] = self.confirmed['files']['spliceai_snv']['path']
        self.assert_rejected('Unexpected SpliceAI assembly/score type')

    def test_resource_filename_for_other_assembly_is_rejected(self):
        self.confirmed['files']['spliceai_indel']['path'] = self.confirmed['files']['spliceai_indel']['path'].replace('hg19', 'hg38')
        self.assert_rejected('Unexpected SpliceAI assembly/score type')

    def test_partial_file_is_rejected(self):
        self.confirmed['files']['spliceai_indel']['path'] += '.partial'
        self.assert_rejected('incomplete file cannot be imported')

    def test_acquisition_size_mismatch_is_rejected(self):
        self.confirmed['files']['spliceai_indel']['bytes'] += 1
        self.assert_rejected('Acquisition byte-size mismatch')

    def test_acquisition_checksum_mismatch_is_rejected(self):
        self.confirmed['files']['spliceai_indel']['sha256'] = '0' * 64
        self.assert_rejected('Acquisition checksum mismatch')

    def test_missing_selected_index_is_rejected(self):
        Path(self.confirmed['files']['spliceai_indel']['path'] + '.tbi').unlink()
        self.assert_rejected('Missing matching index')

    def test_index_symlink_cannot_escape_its_resource_root(self):
        index = Path(self.confirmed['files']['spliceai_indel']['path'] + '.tbi')
        outside = self.folder / 'outside.tbi'
        index.rename(outside)
        index.symlink_to(outside)
        self.assert_rejected('Unexpected shared-cache index root')

    def test_existing_subset_receipt_is_preserved(self):
        receipt = self.folder / 'licensed-caches-hg19-spliceai.ready.json'
        receipt.write_text('original receipt')
        with self.assertRaisesRegex(RuntimeError, 'receipt already exists'):
            self.call_import()
        self.assertEqual(receipt.read_text(), 'original receipt')
        self.assertFalse(self.config.exists())

    def assert_checksum_change_rejected(self, target):
        # Identical bytes with a changed write timestamp still invalidate the verified input.
        sys.path.insert(0, str(ROOT / 'scripts'))
        import deploy_cluster_resources
        real_checksum = deploy_cluster_resources.checksum
        def checksum_then_mutate(path, *args, **kwargs):
            value = real_checksum(path, *args, **kwargs)
            if path == target:
                stat = path.stat()
                os.utime(path, ns=(stat.st_atime_ns, stat.st_mtime_ns + 1_000_000_000))
            return value
        with mock.patch.object(deploy_cluster_resources, 'checksum', side_effect=checksum_then_mutate), \
                mock.patch('subprocess.run', wraps=subprocess.run) as run:
            with self.assertRaisesRegex(RuntimeError, 'Original shared input changed'):
                self.call_import()
        self.assertFalse(any(call.args[0][0] == 'bash' for call in run.call_args_list))
        self.assertFalse(self.config.exists())
        self.assertEqual(list(self.folder.glob('*.ready.json')), [])

    def test_source_change_during_checksum_stops_before_installer(self):
        self.assert_checksum_change_rejected(Path(self.confirmed['files']['spliceai_indel']['path']))

    def test_index_change_during_checksum_stops_before_installer(self):
        self.assert_checksum_change_rejected(Path(self.confirmed['files']['spliceai_indel']['path'] + '.tbi'))

    def test_source_change_during_installer_prevents_receipt(self):
        real_run = subprocess.run
        target = Path(self.confirmed['files']['spliceai_indel']['path'])
        def run_and_tamper(command, **kwargs):
            result = real_run(command, **kwargs)
            if command[0] == 'bash':
                stat = target.stat()
                os.utime(target, ns=(stat.st_atime_ns, stat.st_mtime_ns + 1_000_000_000))
            return result
        with mock.patch('subprocess.run', side_effect=run_and_tamper):
            with self.assertRaisesRegex(RuntimeError, 'Original shared input changed'):
                self.call_import()
        self.assertTrue(self.config.exists())  # Retain this unaccepted config for diagnosis.
        self.assertEqual(list(self.folder.glob('*.ready.json')), [])

    def test_unselected_config_mutation_prevents_completion_receipt(self):
        real_run = subprocess.run
        def run_and_tamper(command, **kwargs):
            result = real_run(command, **kwargs)
            if command[0] == 'bash':
                settings = yaml.safe_load(self.config.read_text())
                settings['primateai_prescore'] = 'unexpected change'
                self.config.write_text(yaml.safe_dump(settings))
            return result
        with mock.patch('subprocess.run', side_effect=run_and_tamper):
            with self.assertRaisesRegex(RuntimeError, 'modified unselected configuration keys'):
                self.call_import()
        self.assertEqual(list(self.folder.glob('*.ready.json')), [])


if __name__ == '__main__':
    unittest.main()
