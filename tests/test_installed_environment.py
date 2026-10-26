"""Installed-package export and child-only CADD environment routing checks."""
import importlib.util
import json
import os
from pathlib import Path
import subprocess
import tempfile
import unittest
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
        result = subprocess.CompletedProcess([], 0, '{"environments":{}}', '')
        with patch.object(EXPORT.subprocess, 'run', return_value=result), \
                patch.object(EXPORT, 'export_environment', side_effect=AssertionError('no export before core exists')):
            with self.assertRaisesRegex(ValueError, 'not installed'):
                EXPORT.export_installation(config_file=self.config, manager='micromamba')

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


if __name__ == '__main__':
    unittest.main()
