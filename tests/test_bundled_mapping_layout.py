"""Mapping storage, explicit input identity and external acquisition boundaries."""
import ast
import hashlib
import importlib.util
import json
import os
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest
from unittest import mock

ROOT = Path(__file__).resolve().parents[1]


def load(name, filename):
    spec = importlib.util.spec_from_file_location(name, ROOT / 'scripts' / filename)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


P = load('layout_prerequisites', 'installer_prerequisites.py')
C = load('layout_config_paths', 'installer_config_paths.py')
B = P.BIOMART
HGNC = 'hgnc_non_alt_loci_set_sha256_890b18743cc8.tsv'


class BundledMappingLayout(unittest.TestCase):
    def setUp(self):
        temp = tempfile.TemporaryDirectory(dir=os.environ.get('PRIVA_TEST_TMPDIR', os.environ.get('TMPDIR')))
        self.addCleanup(temp.cleanup)
        self.folder = Path(temp.name)
        self.checkout = self.folder / 'checkout'
        self.checkout.mkdir()
        self.cache = self.folder / 'cache'
        self.config = {'assembly': 'hg38', 'vep_installed_version': '116',
                       'cache_root': str(self.cache), 'tmp_dir': str(self.folder / 'tmp')}

    def test_every_automatic_mapping_stays_in_repo_data_when_missing(self):
        for assembly, release in (('hg19', '116'), ('GRCh37', '113'), ('hg38', '113'), ('GRCh38', '116')):
            with self.subTest(assembly=assembly, release=release):
                config = dict(self.config, assembly=assembly, vep_installed_version=release)
                selected = B.source_for_assembly(assembly, release)
                target = self.checkout / selected['relative_path']
                before = dict(config)
                self.assertEqual(P.nm_to_enst_path(config, self.checkout), target)
                self.assertEqual(P.mapping_fetch_source(config, self.checkout), selected['source_id'])
                updates = C.portable_updates(config, self.checkout)
                self.assertNotIn('nm_to_enst_tsv', updates)
                self.assertEqual(config, before)
        self.assertFalse(self.cache.exists())
        self.assertFalse((self.checkout / 'data').exists())

    def test_explicit_cache_mappings_and_sidecars_remain_the_original_inputs(self):
        for key in ('nm_to_enst_tsv', 'hgnc_table'):
            with self.subTest(key=key):
                selected = self.cache / 'selected' / key / 'original.tsv'
                selected.parent.mkdir(parents=True)
                selected.write_bytes(b'exact explicitly selected scientific input')
                sidecar = Path(str(selected) + '.provenance.json')
                sidecar.write_bytes(b'{"selected":"original"}\n')
                before = (selected.read_bytes(), selected.stat().st_mtime_ns,
                          sidecar.read_bytes(), sidecar.stat().st_mtime_ns)
                config = dict(self.config, **{key: str(selected)})
                updates = C.portable_updates(config, self.checkout)
                self.assertEqual(dict(config, **updates)[key], str(selected))
                self.assertFalse(any(row['key'] == key for row in C.adoption_plan(config, updates, self.checkout)))
                self.assertEqual(before, (selected.read_bytes(), selected.stat().st_mtime_ns,
                                          sidecar.read_bytes(), sidecar.stat().st_mtime_ns))
                self.assertFalse((self.cache / 'data/hgnc').exists())

    def test_relative_mapping_selection_is_resolved_from_checkout_without_copying(self):
        selected = self.checkout / 'data/hgnc/custom.selected.tsv'
        selected.parent.mkdir(parents=True)
        selected.write_bytes(b'explicit original mapping')
        config = dict(self.config, nm_to_enst_tsv='data/hgnc/custom.selected.tsv')
        updates = C.portable_updates(config, self.checkout)
        self.assertEqual(updates['nm_to_enst_tsv'], str(selected))
        self.assertEqual(C.adoption_plan(config, updates, self.checkout), [])
        self.assertFalse(self.cache.exists())

    def test_adding_ftp_bundle_does_not_stale_explicit_mapping_derived_receipts(self):
        selected = self.cache / 'previous-profile/selected.tsv'
        selected.parent.mkdir(parents=True)
        selected.write_text('\t'.join(B.HEADERS) + '\nENSG1\tENST1\tG\tP\tNM_1\n')
        artifact = self.cache / 'completed-motifs.pkl'
        artifact.write_bytes(b'completed prior derived fixture')
        receipt = Path(str(artifact) + '.nm_to_enst.json')
        receipt.write_text(json.dumps({'mapping_identity': B.mapping_identity(selected),
                                       'artifact_bytes': artifact.stat().st_size,
                                       'artifact_sha256': B.sha256(artifact)}))
        before = (artifact.read_bytes(), artifact.stat().st_mtime_ns, receipt.read_bytes(), receipt.stat().st_mtime_ns)
        bundle = self.checkout / B.FTP.RELATIVE_PATH
        bundle.parent.mkdir(parents=True)
        bundle.write_bytes(b'newly available distinct canonical bundle')
        config = dict(self.config, nm_to_enst_tsv=str(selected))
        resolved = dict(config, **C.portable_updates(config, self.checkout))
        self.assertEqual(P.validate_mapping(resolved, self.checkout), selected)
        B.check_receipt(artifact, selected)
        self.assertEqual(B.artifact_route(artifact=artifact, mapping=selected, cache_root=self.cache,
                                         output_dir=self.cache / 'AlphaMissense', assembly='hg38', kind='motifs'), artifact)
        self.assertEqual(before, (artifact.read_bytes(), artifact.stat().st_mtime_ns,
                                  receipt.read_bytes(), receipt.stat().st_mtime_ns))
        self.assertFalse((self.cache / 'data/hgnc').exists())

    def test_missing_explicit_mapping_fails_before_adoption_or_config_writes(self):
        config = dict(self.config, nm_to_enst_tsv=str(self.cache / 'missing.tsv'))
        updates = C.portable_updates(config, self.checkout)
        with self.assertRaisesRegex(ValueError, 'refusing to substitute'):
            C.adoption_plan(config, updates, self.checkout)
        with self.assertRaisesRegex(ValueError, 'restore that selected input'):
            P.mapping_fetch_source(config, self.checkout)
        self.assertFalse(self.cache.exists())

    def test_workspace_is_external_and_never_a_source_data_destination(self):
        expected = self.cache / 'BioMart' / B.FTP.SOURCE_ID
        self.assertEqual(P.mapping_workspace(self.config, self.checkout), expected)
        temporary = dict(self.config, cache_root=None)
        self.assertEqual(P.mapping_workspace(temporary, self.checkout), self.folder / 'tmp/BioMart' / B.FTP.SOURCE_ID)
        for cache in (str(self.checkout / 'cache'), '.', 'null', None):
            config = dict(self.config, cache_root=cache, tmp_dir='')
            with self.assertRaises(ValueError):
                P.mapping_workspace(config, self.checkout)
        alias = self.folder / 'alias'
        alias.symlink_to(self.checkout, target_is_directory=True)
        with self.assertRaisesRegex(ValueError, 'outside the checkout'):
            P.mapping_workspace(dict(self.config, cache_root=str(alias / 'work')), self.checkout)
        self.assertFalse(expected.exists())

    def test_new_ftp_acquisition_requires_an_explicit_external_workspace(self):
        output = self.folder / B.FTP.BASENAME
        with mock.patch.object(B.FTP, 'install', side_effect=AssertionError('acquisition must not start')):
            with self.assertRaisesRegex(ValueError, '--resource-dir'):
                B.acquire(output, source_id=B.FTP.SOURCE_ID)
            with self.assertRaisesRegex(ValueError, 'outside the checkout'):
                B.acquire(output, source_id=B.FTP.SOURCE_ID, resource_dir=ROOT / 'data/hgnc/raw')
        with mock.patch.object(B.FTP, 'install', return_value={'fixture': True}) as install:
            B.acquire(output, source_id=B.FTP.SOURCE_ID, resource_dir=self.folder / 'work')
            self.assertEqual(install.call_args.args[1], self.folder / 'work')
        with mock.patch.object(B, '_acquire_locked', side_effect=AssertionError('acquisition must not start')):
            with self.assertRaisesRegex(ValueError, 'outside the checkout'):
                B.acquire(output, ROOT / 'data/hgnc/raw.tsv')
        self.assertFalse(output.exists())

    def test_bundled_filename_symlinks_even_within_data_require_full_validation(self):
        selected = self.checkout / B.FTP.RELATIVE_PATH
        selected.parent.mkdir(parents=True)
        payload = selected.with_name('arbitrary.tsv')
        payload.write_bytes(b'not admitted under canonical name')
        selected.symlink_to(payload.name)
        with mock.patch.object(B, 'verify_pinned', side_effect=ValueError('strict external check')) as verify:
            with self.assertRaisesRegex(ValueError, 'strict external'):
                P.validate_mapping(self.config, self.checkout)
            verify.assert_called_once()

    def test_portable_ftp_metadata_stays_strict_for_external_sources(self):
        output = self.folder / B.FTP.BASENAME
        output.write_text('\t'.join(B.HEADERS) + '\nENSG1\tENST1\tG\tP\tNM_1\n')
        record = dict(B.FTP.specification(), metadata_kind='priva_bundled_mapping_snapshot',
                      metadata_schema_version=1, bytes=output.stat().st_size,
                      sha256=hashlib.sha256(output.read_bytes()).hexdigest(), statistics=B.validate_table(output))
        with mock.patch.multiple(B.FTP, EXPECTED_BYTES=record['bytes'], EXPECTED_SHA256=record['sha256'],
                                 EXPECTED_STATISTICS=record['statistics']):
            self.assertEqual(B.FTP.verify_record(output, record, B.validate_table, 'hg38'), record)
            with self.assertRaisesRegex(ValueError, 'metadata schema'):
                B.FTP.verify_record(output, dict(record, metadata_schema_version=True), B.validate_table)
            with self.assertRaisesRegex(ValueError, 'provenance differs'):
                B.FTP.verify_record(output, dict(record, order_policy='sorted'), B.validate_table)
            output.write_text(output.read_text().replace('NM_1', 'NM_2'))
            with self.assertRaisesRegex(ValueError, 'SHA256'):
                B.FTP.verify_record(output, record, B.validate_table)

    def test_hgnc_canonical_and_legacy_defaults_reuse_opaque_bytes_only(self):
        canonical = self.checkout / 'data/hgnc' / HGNC
        canonical.parent.mkdir(parents=True)
        canonical.write_bytes(b'opaque HGNC fixture; do not parse/hash as an installation check')
        legacy = canonical.with_name('non_alt_loci_set.tsv')
        legacy.symlink_to(canonical.name)
        self.assertEqual(C.portable_updates(self.config, self.checkout)['hgnc_table'], str(canonical))
        with mock.patch.object(Path, 'open', side_effect=AssertionError('no bundle content validation')):
            for selected in (str(canonical), 'data/hgnc/' + HGNC, str(legacy), 'data/hgnc/non_alt_loci_set.tsv'):
                self.assertEqual(C.bundled_resource('hgnc_table', selected, self.checkout),
                                 legacy if Path(selected).name == legacy.name else canonical)
        config = dict(self.config, hgnc_table='data/hgnc/non_alt_loci_set.tsv')
        updates = C.portable_updates(config, self.checkout)
        self.assertEqual(updates['hgnc_table'], str(legacy))
        self.assertEqual(C.adoption_plan(config, updates, self.checkout), [])
        self.assertFalse(self.cache.exists())

    def test_registry_and_actual_hgnc_defaults_agree_without_scientific_edits(self):
        registry = json.loads((ROOT / 'cache_resources.json').read_text())['resources']
        self.assertEqual(registry['nm_to_enst']['bundled_paths'],
                         ['data/hgnc/' + name for name in B.PINNED_BASENAMES])
        self.assertIn('data/hgnc/' + HGNC, registry['curated_gene_tables']['bundled_paths'])
        for filename, constant in (('gene_mechanism_common.py', 'DEFAULT_HGNC_TABLE'),
                                   ('build_gene_nonlof_mechanism_cache.py', 'DEFAULT_HGNC_TABLE'),
                                   ('audit_gofcards_source_coverage.py', 'DEFAULT_HGNC')):
            source = (ROOT / 'scripts' / filename).read_text()
            assignment = next(node for node in ast.parse(source).body
                              if isinstance(node, ast.Assign) and any(isinstance(target, ast.Name) and target.id == constant for target in node.targets))
            self.assertEqual(ast.literal_eval(assignment.value.right), HGNC)

    def test_both_missing_source_shell_routes_use_repo_target_and_external_workspace(self):
        for release in ('113', '116'):
            with self.subTest(release=release):
                config = dict(self.config, vep_installed_version=release)
                path = self.folder / ('config-' + release + '.json')
                path.write_text(json.dumps(config))
                events = self.folder / ('events-' + release + '.json')
                # Real selection/workspace CLI; only network acquisition and
                # its immediate full verification are replaced by a fixture.
                script = r'''
source "$1" >/dev/null || exit
BASE_DIR=$TEST_SOURCE
python() {
    if [[ $1 == */biomart_nm_to_enst.py ]]; then
        "$TEST_PYTHON" -c 'import json,sys; from pathlib import Path; a=sys.argv[2:]; p=Path(sys.argv[1]); old=json.loads(p.read_text()) if p.exists() else []; old.append(a); p.write_text(json.dumps(old))' "$TEST_EVENTS" "$@"
        [[ $2 != fetch ]] || return 73
    else "$TEST_PYTHON" "$@"; fi
}
nm_to_enst_mapping_install "$TEST_CONFIG"
'''
                result = subprocess.run(['bash', '--noprofile', '--norc', '-c', script,
                                         'fixture', str(ROOT / 'scripts/install_utils.sh')],
                                        env=dict(os.environ, TEST_PYTHON=sys.executable, TEST_SOURCE=str(self.checkout),
                                                 TEST_CONFIG=str(path), TEST_EVENTS=str(events)),
                                        text=True, capture_output=True, timeout=30)
                self.assertNotEqual(result.returncode, 0)
                argv, = json.loads(events.read_text())
                selected = B.source_for_assembly('hg38', release)
                self.assertEqual(argv[argv.index('--output') + 1], str(self.checkout / selected['relative_path']))
                workspace = self.cache / 'BioMart' / selected['source_id']
                if release == '116':
                    self.assertEqual(argv[argv.index('--resource-dir') + 1], str(workspace))
                    self.assertNotIn('--raw-output', argv)
                else:
                    self.assertTrue(Path(argv[argv.index('--raw-output') + 1]).is_relative_to(workspace))
                self.assertEqual(json.loads(path.read_text()), config)
                self.assertFalse((self.cache / 'data/hgnc').exists())


if __name__ == '__main__':
    unittest.main()
