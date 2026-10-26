"""Release selection and immutable mapping/cache boundaries for new installs."""
import hashlib
import importlib.util
import json
import os
from pathlib import Path
import tempfile
import unittest
from unittest import mock

ROOT = Path(__file__).resolve().parents[1]
spec = importlib.util.spec_from_file_location('mapping_profile_prerequisites', ROOT / 'scripts/installer_prerequisites.py')
P = importlib.util.module_from_spec(spec)
spec.loader.exec_module(P)
B = P.BIOMART


class MappingReleaseSelection(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory(dir=os.environ.get('PRIVA_TEST_TMPDIR', os.environ.get('TMPDIR')))
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)
        self.env = mock.patch.dict(os.environ, {'PRIVA_INSTALL_INSPECT_ROOTS': str(self.root)})
        self.env.start()
        self.addCleanup(self.env.stop)

    def test_exact_116_selects_complete_ftp_source_for_grch38_aliases(self):
        for assembly in ('hg38', 'GRCh38', 'GRCh38.p14'):
            for release in ('116', 116):
                source = B.source_for_assembly(assembly, release)
                self.assertEqual(source['source_id'], B.FTP.SOURCE_ID)
                self.assertEqual(source['expected_bytes'], 28785905)
                self.assertEqual(source['expected_sha256'], B.FTP.EXPECTED_SHA256)

    def test_historical113_and_no_release_keep_original_source(self):
        for release in (None, '113', 113):
            self.assertEqual(B.source_for_assembly('hg38', release)['source_id'], B.GRCH38_SOURCE_ID)

    def test_grch37_mart116_annotation75_is_separate(self):
        for release in (None, '113', '116'):
            for assembly in ('hg19', 'GRCh37', 'GRCh37.p13'):
                source = B.source_for_assembly(assembly, release)
                self.assertEqual(source['source_id'], B.SOURCE_ID)
                self.assertEqual((source['mart_release'], source['annotation_base_release']), (116, 75))

    def test_unsupported_major_has_no_automatic_latest_fallback(self):
        for release in ('114', '115', '117', 'latest', '116.3', True):
            with self.assertRaisesRegex(ValueError, 'explicit.*nm_to_enst_tsv'):
                B.source_for_assembly('hg38', release)

    def test_null_116_default_uses_repository_data_even_before_acquisition(self):
        config = {'assembly': 'hg38', 'vep_installed_version': '116', 'cache_root': str(self.root / 'payloads'), 'nm_to_enst_tsv': None}
        before = dict(config)
        target = self.root / 'data/hgnc' / B.FTP.BASENAME
        self.assertEqual(P.nm_to_enst_path(config, self.root), target)
        self.assertEqual(P.mapping_fetch_source(config, self.root), B.FTP.SOURCE_ID)
        self.assertEqual(config, before)
        self.assertFalse(target.parent.exists())
        self.assertFalse((self.root / 'payloads').exists())
        with mock.patch.object(P, 'validate_selected_clinvar_stats'):
            P.validate_main(config, self.root, allow_missing_bundled_mapping=True)

    def test_genuine_existing_bundle_reuses_only_nonzero_existence(self):
        for assembly, release in (('hg19', '116'), ('hg38', '113'), ('hg38', '116')):
            selected = self.root / B.source_for_assembly(assembly, release)['relative_path']
            selected.parent.mkdir(parents=True, exist_ok=True)
            selected.write_bytes(b'opaque bundled bytes; intentionally not schema-validated')
            config = {'assembly': assembly, 'vep_installed_version': release, 'cache_root': str(self.root / 'cache')}
            with mock.patch.object(B, 'verify_pinned', side_effect=AssertionError('extra validation forbidden')):
                self.assertEqual(P.validate_mapping(config, self.root), selected)
            self.assertFalse((self.root / 'cache').exists())
            selected.write_bytes(b'')
            with self.assertRaises(ValueError):
                P.validate_mapping(config, self.root)

    def test_external_symlink_and_external_ftp_resource_keep_full_validation(self):
        selected = self.root / B.GRCH38_RELATIVE_PATH
        selected.parent.mkdir(parents=True)
        external = self.root / 'external.tsv'
        external.write_bytes(b'external downloaded bytes')
        selected.symlink_to(external)
        config = {'assembly': 'hg38', 'nm_to_enst_tsv': str(selected)}
        with mock.patch.object(B, 'verify_pinned', side_effect=ValueError('external must validate')) as strict:
            with self.assertRaisesRegex(ValueError, 'external must validate'):
                P.validate_mapping(config, self.root)
            strict.assert_called_once()
        ftp = self.root / B.FTP.BASENAME
        ftp.write_bytes(b'explicit external FTP snapshot, outside canonical bundle path')
        config['nm_to_enst_tsv'] = str(ftp)
        with mock.patch.object(B, 'verify_pinned', side_effect=ValueError('FTP must validate')):
            with self.assertRaisesRegex(ValueError, 'FTP must validate'):
                P.validate_mapping(config, self.root)

    def test_historical_source_root_default_and_relative_cache_root(self):
        config = {'assembly': 'hg38', 'vep_installed_version': '113'}
        self.assertEqual(P.nm_to_enst_path(config, self.root), self.root / B.GRCH38_RELATIVE_PATH)
        self.assertEqual(P.mapping_fetch_source(config, self.root), B.GRCH38_SOURCE_ID)
        config['cache_root'] = 'payloads'
        self.assertEqual(P.nm_to_enst_path(config, self.root), self.root / B.GRCH38_RELATIVE_PATH)

    def test_explicit_path_bypasses_release_selection_but_not_assembly(self):
        table = self.root / 'selected.tsv'
        table.write_text('\t'.join(B.HEADERS) + '\nENSG1\tENST1\tG\tP\tNM_1\n')
        config = {'assembly': 'hg38', 'vep_installed_version': '999', 'nm_to_enst_tsv': str(table)}
        self.assertEqual(P.validate_mapping(config, self.root), table)
        config.update(assembly='hg19', nm_to_enst_tsv=str(self.root / B.FTP.BASENAME))
        with self.assertRaisesRegex(ValueError, 'not compatible'):
            P.nm_to_enst_path(config, self.root)

    def test_missing_explicit_never_authorizes_fetch_even_at_default_destination(self):
        config = {'assembly': 'hg38', 'vep_installed_version': '116', 'cache_root': str(self.root)}
        target = P.nm_to_enst_path(config, self.root)
        for release in ('116', '999'):
            config.update(vep_installed_version=release, nm_to_enst_tsv=str(target))
            with self.assertRaisesRegex(ValueError, 'Explicit nm_to_enst_tsv is missing'):
                P.mapping_fetch_source(config, self.root)
        self.assertFalse(target.exists())

    def test_order_and_multiplicity_change_reject_existing_derived_receipt(self):
        header = '\t'.join(B.HEADERS) + '\n'
        a, b = 'ENSG1\tENST1\tG\tP\tNM_1\n', 'ENSG1\tENST2\tG\tP\tNM_1\n'
        prior, updated = self.root / 'prior.tsv', self.root / 'updated.tsv'
        prior.write_text(header + a + b)
        artifact = self.root / 'derived.pkl'
        artifact.write_bytes(b'prior scientific output remains unchanged')
        receipt = {'mapping_identity': B.mapping_identity(prior), 'artifact_bytes': artifact.stat().st_size, 'artifact_sha256': hashlib.sha256(artifact.read_bytes()).hexdigest()}
        Path(str(artifact) + '.nm_to_enst.json').write_text(json.dumps(receipt))
        B.check_receipt(artifact, prior)
        for rows in (b + a, a + b + a):
            updated.write_text(header + rows)
            with self.assertRaisesRegex(ValueError, 'Cannot reuse or overwrite'):
                B.check_receipt(artifact, updated)
        self.assertEqual(artifact.read_bytes(), b'prior scientific output remains unchanged')

    def test_validator_retains_blank_nm_and_duplicate_rows_in_original_order(self):
        table, output = self.root / 'raw.tsv', self.root / 'retained.tsv'
        a = 'ENSG1\tENST1\tG\tP\tNM_1\n'
        body = '\t'.join(B.HEADERS) + '\n' + a + 'ENSG2\tENST2\t\t\t\n' + a
        table.write_text(body + '[success]\n')
        counts = B.validate_table(table, raw=True, output=output)
        self.assertEqual(output.read_text(), body)
        self.assertEqual((counts['rows'], counts['rows_without_nm'], counts['duplicate_full_rows'], counts['duplicate_nm_enst_rows']), (3, 1, 1, 1))


if __name__ == '__main__':
    unittest.main()
