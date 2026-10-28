"""Installed-package export and child-only CADD environment routing checks."""
import importlib.util
import hashlib
import json
import os
from pathlib import Path
import subprocess
import tempfile
import unittest
import uuid
from unittest.mock import patch

import yaml

ROOT = Path(__file__).resolve().parents[1]
SPEC = importlib.util.spec_from_file_location('installed_environment_test', ROOT / 'scripts/installed_environment.py')
EXPORT = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(EXPORT)
RUNTIME_SPEC = importlib.util.spec_from_file_location('cadd_environment_route_test', ROOT / 'scripts/cadd_runtime.py')
RUNTIME = importlib.util.module_from_spec(RUNTIME_SPEC)
RUNTIME_SPEC.loader.exec_module(RUNTIME)


class InstalledEnvironment(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory(dir=os.environ.get('PRIVA_TEST_TMPDIR', os.environ.get('TMPDIR')))
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)
        self.config = self.root / 'config.yml'
        self.config.write_text(yaml.safe_dump(dict(cache_root=str(self.root))))
        self.export = dict(name='priva_acmg', channels=['conda-forge', 'bioconda'],
            dependencies=['python=3.12.3=build', 'numpy=1.26.4=build', {'pip': ['actual-package==1.2.3']}],
            prefix='/host/envs/priva_acmg', variables={'PRIVATE_SETTING': 'must not be exported'})

    def invoke(self, document=None, code=0):
        result = subprocess.CompletedProcess([], code, yaml.safe_dump(self.export if document is None else document), 'failed export')
        with patch.object(EXPORT.subprocess, 'run', return_value=result) as call:
            output = EXPORT.export_environment(config_file=self.config, manager='/selected/micromamba')
        self.assertEqual(call.call_args.args[0], ['/selected/micromamba', 'env', 'export', '--name', 'priva_acmg'])
        return output

    def test_export_uses_installed_versions_and_omits_host_settings(self):
        output = self.invoke()
        document = yaml.safe_load(output.read_text())
        self.assertEqual(set(document), {'name', 'channels', 'dependencies'})
        self.assertEqual(document['dependencies'], self.export['dependencies'])
        before = output.read_bytes(), output.stat().st_mtime_ns
        self.invoke()
        self.assertEqual(before, (output.read_bytes(), output.stat().st_mtime_ns))

    def test_failed_export_preserves_previous_snapshot(self):
        output = self.invoke()
        before = output.read_bytes()
        with self.assertRaisesRegex(ValueError, 'could not export'):
            self.invoke(code=17)
        self.assertEqual(output.read_bytes(), before)

    def test_different_environment_or_empty_export_cannot_publish(self):
        for document in (dict(self.export, name='another_env'), dict(self.export, dependencies=[])):
            with self.assertRaises(ValueError):
                self.invoke(document)
        self.assertFalse((self.root / 'priva_acmg.installed.yml').exists())

    def test_link_or_unrelated_destination_is_retained(self):
        output = self.root / 'priva_acmg.installed.yml'
        external = self.root / 'other.yml'
        external.write_text('original unrelated content')
        output.symlink_to(external)
        with self.assertRaisesRegex(ValueError, 'linked'):
            self.invoke()
        self.assertEqual(external.read_text(), 'original unrelated content')
        output.unlink()
        output.write_text('unrelated: content\n')
        with self.assertRaisesRegex(ValueError, 'retained'):
            self.invoke()
        self.assertEqual(output.read_text(), 'unrelated: content\n')

    def test_only_environment_variables_are_excluded_not_pip_packages(self):
        document = EXPORT.package_document(self.export, 'priva_acmg')
        self.assertEqual(document['dependencies'][-1], {'pip': ['actual-package==1.2.3']})
        self.assertNotIn('variables', document)

    def test_explicit_cadd_prefix_exports_packages_without_activation(self):
        prefix = self.root / 'cadd-python2'
        (prefix / 'conda-meta').mkdir(parents=True)
        (prefix / 'conda-meta/history').write_text('completed fixture')
        document = dict(name=prefix.name, prefix=str(prefix), channels=['conda-forge'], dependencies=['python=2.7.13=0'])
        result = subprocess.CompletedProcess([], 0, yaml.safe_dump(document), '')
        with patch.object(EXPORT.subprocess, 'run', return_value=result) as call:
            output = EXPORT.export_environment(config_file=self.config, manager='/selected/micromamba',
                                                name=prefix.name, prefix=prefix)
        self.assertEqual(call.call_args.args[0], ['/selected/micromamba', 'env', 'export', '--prefix', str(prefix)])
        self.assertEqual(yaml.safe_load(output.read_text()), {key: value for key, value in document.items() if key != 'prefix'})

    def test_micromamba_unnamed_prefix_requires_exact_reported_prefix(self):
        prefix = self.root / 'cadd-python2'
        document = dict(name=None, prefix=str(prefix), channels=['bioconda'], dependencies=['python=2.7.13=0'])
        self.assertEqual(EXPORT.package_document(document, prefix.name, selected_prefix=prefix)['name'], prefix.name)
        with self.assertRaisesRegex(ValueError, 'different environment prefix'):
            EXPORT.package_document(dict(document, prefix=str(self.root / 'other')), prefix.name, selected_prefix=prefix)
        with self.assertRaisesRegex(ValueError, 'different environment prefix'):
            EXPORT.package_document(dict(document, name=prefix.name, prefix=str(self.root / 'other')),
                                    prefix.name, selected_prefix=prefix)
        with self.assertRaisesRegex(ValueError, 'different environment'):
            EXPORT.package_document(document, prefix.name)

    def test_installation_exports_only_completed_cadd_prefixes(self):
        environments = {}
        for name, present in (('cadd-python2', True), ('cadd-esm-python', True), ('cadd-mmsplice-python', False)):
            prefix = self.root / name
            if present:
                (prefix / 'conda-meta').mkdir(parents=True)
                (prefix / 'conda-meta/history').write_text('fixture installed packages')
            environments[name] = {'prefix': str(prefix)}
        result = subprocess.CompletedProcess([], 0, json.dumps({'environments': environments}), '')
        with patch.object(EXPORT.subprocess, 'run', return_value=result), \
                patch.object(EXPORT, 'export_environment', return_value=self.root / 'export.yml') as export:
            self.assertEqual(len(EXPORT.export_installation(config_file=self.config, manager='micromamba')), 3)
        self.assertEqual([call.kwargs.get('name', 'priva_acmg') for call in export.call_args_list],
                         ['priva_acmg', 'cadd-python2', 'cadd-esm-python'])

    def test_installation_does_not_claim_an_absent_cadd_environment(self):
        self.config.write_text(yaml.safe_dump(dict(cache_root=str(self.root), cadd_profile='/selected/profile.json')))
        result = subprocess.CompletedProcess([], 0, '{"environments":{}}', '')
        with patch.object(EXPORT.subprocess, 'run', return_value=result), \
                patch.object(EXPORT, 'export_environment', side_effect=AssertionError('no export before core exists')):
            with self.assertRaisesRegex(ValueError, 'not installed'):
                EXPORT.export_installation(config_file=self.config, manager='micromamba')

    def test_installation_without_profile_exports_priva_when_no_cadd_runtime(self):
        # Upstream CADD.sh needs no separate scorer, so its absence is not an error.
        for profile in (None, '', 'null'):
            with self.subTest(profile=profile):
                self.config.write_text(yaml.safe_dump(dict(cache_root=str(self.root), cadd_profile=profile)))
                result = subprocess.CompletedProcess([], 0, '{"environments":{"cadd-python2":{"prefix":"/absent"}}}', '')
                with patch.object(EXPORT.subprocess, 'run', return_value=result), \
                        patch.object(EXPORT, 'export_environment', return_value=self.root / 'priva_acmg.installed.yml') as export:
                    self.assertEqual(len(EXPORT.export_installation(config_file=self.config, manager='micromamba')), 1)
                self.assertEqual(export.call_count, 1)
                self.assertNotIn('prefix', export.call_args.kwargs)

    def test_explicit_verified_container_exports_no_fictitious_mamba_environment(self):
        result = subprocess.CompletedProcess([], 0, '{"backend":"container","environments":{}}', '')
        with patch.object(EXPORT.subprocess, 'run', return_value=result), \
                patch.object(EXPORT, 'export_environment', return_value=self.root / 'priva_acmg.installed.yml') as export:
            self.assertEqual(len(EXPORT.export_installation(config_file=self.config, manager='micromamba')), 1)
        self.assertEqual(export.call_count, 1)
        self.assertNotIn('prefix', export.call_args.kwargs)

    def test_unknown_runtime_backend_is_not_exported(self):
        result = subprocess.CompletedProcess([], 0, '{"backend":"unknown","environments":{}}', '')
        with patch.object(EXPORT.subprocess, 'run', return_value=result), \
                patch.object(EXPORT, 'export_environment', side_effect=AssertionError('no export for unknown backend')):
            with self.assertRaisesRegex(ValueError, 'Unknown'):
                EXPORT.export_installation(config_file=self.config, manager='micromamba')

    def test_cadd_child_preserves_main_environment_on_success_and_failure(self):
        # A real wrapper and Bash process execute a labelled interpreter fixture;
        # this checks context/exit behavior, not scientific CADD calculation.
        folder = self.root / 'cadd'
        folder.mkdir()
        prefix = folder / 'envs/cadd-python2'
        (prefix / 'bin').mkdir(parents=True)
        python = prefix / 'bin/python'
        python.write_text('#!/usr/bin/env bash\n'
            'printf "CHILD_PREFIX=%s\\n" "$CONDA_PREFIX"\n'
            'printf "CHILD_PYTHONPATH=%s\\n" "${PYTHONPATH-unset}"\n'
            'exit "$1"\n')
        python.chmod(0o755)
        routes = RUNTIME.install_cadd_entry_points(folder, native=True, names=('cadd-python2',))
        program = routes['entry_points']['cadd-python2']
        body = ('before_path=$PATH\n'
                '"$1" "$2"\n'
                'status=$?\n'
                '[[ $CONDA_PREFIX == /fixture/priva_acmg && $PYTHONPATH == main_modules && $PATH == "$before_path" ]] || exit 99\n'
                'printf "MAIN_RESTORED\\n"\n'
                'exit "$status"\n')
        for code in (0, 72):
            with self.subTest(code=code):
                result = subprocess.run(['bash', '--noprofile', '--norc', '-c', body, 'test', program, str(code)],
                    env=dict(os.environ, CONDA_PREFIX='/fixture/priva_acmg', PYTHONPATH='main_modules'),
                    capture_output=True, text=True)
                self.assertEqual(result.returncode, code, result.stderr)
                self.assertIn('CHILD_PREFIX=' + str(prefix), result.stdout)
                self.assertIn('CHILD_PYTHONPATH=unset', result.stdout)
                self.assertIn('MAIN_RESTORED', result.stdout)


class ManagedInstalledEnvironment(unittest.TestCase):
    """Exercise real committed annex pointers without writing through them."""

    def setUp(self):
        self.temp = tempfile.TemporaryDirectory(dir=os.environ.get('PRIVA_TEST_TMPDIR', os.environ.get('TMPDIR')))
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)
        self.repo = self.root / 'dataset'
        self.repo.mkdir()
        self.cache = self.repo / 'caches'
        self.cache.mkdir()
        self.config = self.root / 'config.yml'
        self.config.write_text(yaml.safe_dump(dict(cache_root=str(self.cache))))
        self.export = dict(name='priva_acmg', channels=['conda-forge'],
                           dependencies=['python=3.12.3=build', {'pip': ['actual-package==1.2.3']}],
                           prefix='/host/envs/priva_acmg', variables={'PRIVATE': 'omit'})
        self.payload = yaml.safe_dump(EXPORT.package_document(self.export, 'priva_acmg'),
                                      sort_keys=False).encode('utf-8')
        self.real_run = subprocess.run
        self.git('init', '--quiet')
        self.git('config', 'user.name', 'Installed export fixture')
        self.git('config', 'user.email', 'installed-export@example.invalid')
        self.git('config', 'annex.uuid', str(uuid.uuid4()))
        self.output = self.cache / 'priva_acmg.installed.yml'
        self.target = self.object(self.payload)
        self.output.symlink_to(os.path.relpath(self.target, self.cache))
        self.commit_pointer()

    def git(self, *args, input=None):
        result = self.real_run(['git', '-C', str(self.repo), '-c', 'core.hooksPath=/dev/null',
                               '-c', 'commit.gpgSign=false', *args], input=input,
                              capture_output=True, env=dict(os.environ, GIT_CONFIG_NOSYSTEM='1',
                                                            GIT_CONFIG_GLOBAL='/dev/null'))
        self.assertEqual(result.returncode, 0, result.stderr.decode('utf-8', errors='replace'))
        return result.stdout

    def object(self, payload, backend='MD5E'):
        digest = (hashlib.md5(payload, usedforsecurity=False).hexdigest() if backend == 'MD5E'
                  else hashlib.sha256(payload).hexdigest())
        key = backend + '-s' + str(len(payload)) + '--' + digest + '.yml'
        target = self.repo / '.git/annex/objects/fixture/fixture' / key / key
        target.parent.mkdir(parents=True, exist_ok=True)
        target.write_bytes(payload)
        return target

    def commit_pointer(self):
        self.git('add', '--', str(self.output.relative_to(self.repo)))
        self.git('commit', '--quiet', '-m', 'Record installed package pointer')

    @staticmethod
    def identity(value):
        # Reading a valid export can update atime; protect the actual raw-seven
        # link/object identity and content rather than access-time bookkeeping.
        return (value.st_dev, value.st_ino, value.st_mode, value.st_size,
                value.st_mtime_ns, value.st_ctime_ns, value.st_nlink)

    def protected(self):
        return (self.identity(self.output.lstat()), os.readlink(self.output),
                self.identity(self.target.lstat()), self.target.read_bytes(),
                tuple(sorted(path.name for path in self.cache.iterdir())))

    def invoke(self, *, action=None, document=None, code=0):
        calls = []

        def run(argv, **kwargs):
            if argv[0] == '/selected/micromamba':
                calls.append(argv)
                if action:
                    action()
                return subprocess.CompletedProcess(argv, code,
                    yaml.safe_dump(self.export if document is None else document), 'fixture export failure')
            return self.real_run(argv, **kwargs)

        with patch.object(EXPORT.subprocess, 'run', side_effect=run):
            result = EXPORT.export_environment(config_file=self.config, manager='/selected/micromamba')
        self.assertEqual(calls, [['/selected/micromamba', 'env', 'export', '--name', 'priva_acmg']])
        return result

    def test_identical_committed_md5_export_is_read_only(self):
        before = self.protected()
        self.assertEqual(self.invoke(), self.output)
        self.assertEqual(self.protected(), before)

    def test_identical_committed_sha256_export_is_read_only(self):
        self.target = self.object(self.payload, 'SHA256E')
        self.output.unlink()
        self.output.symlink_to(os.path.relpath(self.target, self.cache))
        self.commit_pointer()
        before = self.protected()
        self.invoke()
        self.assertEqual(self.protected(), before)

    def test_changed_installed_packages_preserve_link_and_object(self):
        before = self.protected()
        with self.assertRaisesRegex(ValueError, 'differs.*retained'):
            self.invoke(document=dict(self.export, dependencies=['python=3.13.0=changed']))
        self.assertEqual(self.protected(), before)

    def test_failed_manager_preserves_committed_export(self):
        before = self.protected()
        with self.assertRaisesRegex(ValueError, 'could not export'):
            self.invoke(code=17)
        self.assertEqual(self.protected(), before)

    def test_invalid_export_schema_is_retained(self):
        self.target = self.object(b'name: priva_acmg\nchannels: []\ndependencies: []\nforeign: true\n')
        self.output.unlink()
        self.output.symlink_to(os.path.relpath(self.target, self.cache))
        self.commit_pointer()
        before = self.protected()
        with self.assertRaisesRegex(ValueError, 'retained'):
            self.invoke()
        self.assertEqual(self.protected(), before)

    def test_uncommitted_literal_link_is_refused(self):
        self.target = self.object(self.payload, 'SHA256E')
        self.output.unlink()
        self.output.symlink_to(os.path.relpath(self.target, self.cache))
        before = self.protected()
        with self.assertRaisesRegex(ValueError, 'linked'):
            self.invoke()
        self.assertEqual(self.protected(), before)

    def test_foreign_repository_object_is_refused(self):
        other = self.root / 'foreign/annex/objects' / self.target.parent.name / self.target.name
        other.parent.mkdir(parents=True)
        other.write_bytes(self.payload)
        self.output.unlink()
        self.output.symlink_to(other)
        self.commit_pointer()
        before = self.identity(self.output.lstat()), other.read_bytes(), self.identity(other.stat())
        with self.assertRaisesRegex(ValueError, 'linked'):
            self.invoke()
        self.assertEqual((self.identity(self.output.lstat()), other.read_bytes(), self.identity(other.stat())), before)

    def test_dangling_committed_link_is_refused(self):
        self.target.unlink()
        before = self.identity(self.output.lstat()), os.readlink(self.output)
        with self.assertRaisesRegex(ValueError, 'linked'):
            self.invoke()
        self.assertEqual((self.identity(self.output.lstat()), os.readlink(self.output)), before)
        self.assertFalse(self.target.exists())

    def test_corrupt_annex_object_is_refused_without_repair(self):
        self.target.write_bytes(self.payload.replace(b'python', b'PYTHON'))
        before = self.protected()
        with self.assertRaisesRegex(ValueError, 'linked'):
            self.invoke()
        self.assertEqual(self.protected(), before)

    def test_annex_object_symlink_is_refused(self):
        other = self.object(self.payload, 'SHA256E')
        self.target.unlink()
        self.target.symlink_to(os.path.relpath(other, self.target.parent))
        before = self.identity(self.target.lstat()), os.readlink(self.target), self.identity(other.stat()), other.read_bytes()
        with self.assertRaisesRegex(ValueError, 'linked'):
            self.invoke()
        self.assertEqual((self.identity(self.target.lstat()), os.readlink(self.target), self.identity(other.stat()), other.read_bytes()), before)

    def test_missing_annex_owner_is_refused(self):
        self.git('config', '--unset', 'annex.uuid')
        before = self.protected()
        with self.assertRaisesRegex(ValueError, 'linked'):
            self.invoke()
        self.assertEqual(self.protected(), before)

    def test_link_retarget_during_manager_export_is_refused(self):
        other = self.object(self.payload, 'SHA256E')
        target_before = self.identity(self.target.stat()), self.target.read_bytes()

        def retarget():
            self.output.unlink()
            self.output.symlink_to(os.path.relpath(other, self.cache))

        with self.assertRaisesRegex(ValueError, 'linked'):
            self.invoke(action=retarget)
        self.assertEqual((self.identity(self.target.stat()), self.target.read_bytes()), target_before)
        self.assertEqual(self.output.resolve(), other)

    def test_object_change_during_manager_export_is_refused(self):
        def change():
            self.target.write_bytes(self.payload.replace(b'python', b'PYTHON'))

        with self.assertRaisesRegex(ValueError, 'linked'):
            self.invoke(action=change)
        self.assertEqual(self.target.read_bytes(), self.payload.replace(b'python', b'PYTHON'))

    def test_owner_change_during_manager_export_is_refused(self):
        before = self.protected()
        with self.assertRaisesRegex(ValueError, 'changed.*retained'):
            self.invoke(action=lambda: self.git('config', 'annex.uuid', str(uuid.uuid4())))
        self.assertEqual(self.protected(), before)

    def test_committed_pointer_change_during_manager_export_is_refused(self):
        before = self.protected()

        def change_pointer():
            oid = self.git('hash-object', '-w', '--stdin', input=b'foreign-export.yml').strip().decode('ascii')
            self.git('update-index', '--cacheinfo', '120000,' + oid + ',' + str(self.output.relative_to(self.repo)))
            self.git('commit', '--quiet', '-m', 'Change only committed export pointer')

        with self.assertRaisesRegex(ValueError, 'linked'):
            self.invoke(action=change_pointer)
        self.assertEqual(self.protected(), before)

    def test_unrelated_commit_with_unchanged_pointer_allows_read_only_reuse(self):
        before = self.protected()
        head = self.git('rev-parse', 'HEAD')

        def other_commit():
            (self.repo / 'metadata.txt').write_text('unrelated completed evidence\n')
            self.git('add', '--', 'metadata.txt')
            self.git('commit', '--quiet', '-m', 'Record unrelated metadata')

        self.invoke(action=other_commit)
        self.assertNotEqual(self.git('rev-parse', 'HEAD'), head)
        self.assertEqual(self.protected(), before)

    def test_object_replacement_during_final_pointer_check_is_refused(self):
        real_git = EXPORT._export_git
        reads = 0

        def checked_git(directory, *args):
            nonlocal reads
            result = real_git(directory, *args)
            if args[:1] == ('ls-tree',):
                reads += 1
                if reads == 2:
                    replacement = self.target.with_name('replacement')
                    replacement.write_bytes(self.payload)
                    os.replace(replacement, self.target)
            return result

        with patch.object(EXPORT, '_export_git', side_effect=checked_git):
            with self.assertRaisesRegex(ValueError, 'linked'):
                self.invoke()
        self.assertEqual(self.target.read_bytes(), self.payload)

    def test_new_foreign_link_during_manager_export_is_refused(self):
        self.output.unlink()
        target_before = self.identity(self.target.stat()), self.target.read_bytes()
        with self.assertRaisesRegex(ValueError, 'linked'):
            self.invoke(action=lambda: self.output.symlink_to(self.target))
        self.assertEqual((self.identity(self.target.stat()), self.target.read_bytes()), target_before)
        self.assertTrue(self.output.is_symlink())


if __name__ == '__main__':
    unittest.main()
