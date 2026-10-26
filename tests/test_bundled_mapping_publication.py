"""Real publication boundaries with tiny mocked publisher exports/projections."""
from contextlib import ExitStack
import hashlib
import importlib.util
import io
import json
import os
from pathlib import Path
import shutil
import tempfile
import unittest
from unittest import mock

from test_biomart_nm_to_enst import BODY, REGISTRY, DATASETS, REGISTRY38, DATASETS38

ROOT = Path(__file__).resolve().parents[1]


class BundledMappingPublication(unittest.TestCase):
    def setUp(self):
        temporary = tempfile.TemporaryDirectory(dir=os.environ.get('PRIVA_TEST_TMPDIR', os.environ.get('TMPDIR')))
        self.addCleanup(temporary.cleanup)
        self.folder = Path(temporary.name)
        self.checkout = self.folder / 'checkout'
        scripts = self.checkout / 'scripts'
        scripts.mkdir(parents=True)
        for name in ('biomart_nm_to_enst.py', 'biomart_ftp_nm_to_enst.py'):
            shutil.copyfile(ROOT / 'scripts' / name, scripts / name)
        spec = importlib.util.spec_from_file_location('publication_mapping', scripts / 'biomart_nm_to_enst.py')
        self.mapping = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(self.mapping)
        self.work = self.folder / 'external-workspace'
        self.work.mkdir()

    def assert_portable(self, record):
        self.assertEqual(record['metadata_kind'], 'priva_bundled_mapping_snapshot')
        for key in ('slurm_job_id', 'producer_script_path', 'retrieved_utc', 'registered_utc',
                    'registration_producer_path', 'provenance_basis', 'operation'):
            self.assertNotIn(key, record)
        for key in ('scientific_qualified', 'production_selected', 'original_biomart_export_identity'):
            self.assertIs(record[key], False)
        for key in ('source_id', 'assembly', 'mart_release', 'annotation_base_release',
                    'header', 'attributes', 'query_sha256', 'bytes', 'sha256', 'order_policy', 'multiplicity_policy'):
            self.assertIn(key, record)
        self.assertEqual(record['bytes'], len(BODY))
        self.assertEqual(record['sha256'], hashlib.sha256(BODY).hexdigest())

    def service(self, assembly, *, canonical=True):
        module = self.mapping
        source = module.source_for_assembly(assembly)
        output = (self.checkout / source['relative_path'] if canonical
                  else self.folder / 'external' / source['basename'])
        raw = self.work / (source['basename'] + '.raw.tsv')
        metadata = module.parse_metadata(REGISTRY if assembly == 'hg19' else REGISTRY38,
                                         DATASETS if assembly == 'hg19' else DATASETS38,
                                         source['source_id'])
        constants = ({'EXPECTED_SNAPSHOT_BYTES': len(BODY), 'EXPECTED_SNAPSHOT_SHA256': hashlib.sha256(BODY).hexdigest()}
                     if assembly == 'hg19' else
                     {'GRCH38_EXPECTED_SNAPSHOT_BYTES': len(BODY), 'GRCH38_EXPECTED_SNAPSHOT_SHA256': hashlib.sha256(BODY).hexdigest()})
        with mock.patch.multiple(module, **constants), \
                mock.patch.object(module, 'source_metadata', return_value=metadata), \
                mock.patch.object(module, '_request', side_effect=lambda *a: io.BytesIO(BODY + b'[success]\n')):
            record = module.acquire(output, raw, source_id=source['source_id'])
            self.assertEqual(module.verify_pinned(output, assembly), record)
        return output, raw, record, constants

    def ftp(self, *, canonical=True):
        module, ftp = self.mapping, self.mapping.FTP
        output = (self.checkout / ftp.RELATIVE_PATH if canonical
                  else self.folder / 'external' / ftp.BASENAME)

        def acquire(raw, name, timeout):
            content = ''.join(f'{bsd[0]} {bsd[1]} {name}\n' for name, (_, bsd, _) in ftp.PINS.items() if bsd)
            (raw / name).write_text(content if name == 'CHECKSUMS' else 'mocked pinned raw input\n')

        def export(database, destination, audit):
            destination.write_bytes(BODY)

        table = self.folder / 'expected.tsv'
        table.write_bytes(BODY)
        constants = {'EXPECTED_BYTES': len(BODY), 'EXPECTED_SHA256': hashlib.sha256(BODY).hexdigest(),
                     'EXPECTED_STATISTICS': module.validate_table(table)}
        with ExitStack() as stack:
            stack.enter_context(mock.patch.multiple(ftp, **constants))
            for name, value in (('acquire', acquire), ('read_schema', lambda *a: {table: [] for table, _ in ftp.PROJECTIONS.values()}),
                                ('metadata', lambda *a: None), ('load_projection', lambda *a: None),
                                ('relationship_audit', lambda *a: {}), ('export_candidate', export)):
                stack.enter_context(mock.patch.object(ftp, name, side_effect=value))
            stack.enter_context(mock.patch.object(ftp, 'urlopen', side_effect=AssertionError('no network in fixtures')))
            record = module.acquire(output, source_id=ftp.SOURCE_ID, resource_dir=self.work)
            self.assertEqual(ftp.verify(output, module.validate_table, 'hg38'), record)
        return output, record, constants

    def test_new_canonical_service_tables_have_portable_metadata_and_external_full_receipts(self):
        for assembly in ('hg19', 'hg38'):
            with self.subTest(assembly=assembly):
                output, raw, record, _ = self.service(assembly)
                self.assert_portable(record)
                self.assertIn('endpoint', record)
                self.assertIn('metadata_before', record)
                self.assertIn('metadata_after', record)
                full, = self.work.glob(output.name + '.acquisition.*.json')
                runtime = json.loads(full.read_text())
                self.assertIn('producer_script_path', runtime)
                self.assertIn('retrieved_utc', runtime)
                self.assertIn('slurm_job_id', runtime)
                self.assertTrue(full.is_relative_to(raw.parent))
                self.assertEqual(output.read_bytes(), BODY)

    def test_new_external_service_acquisition_keeps_its_original_full_sidecar_behavior(self):
        output, _, record, _ = self.service('hg19', canonical=False)
        self.assertIn('producer_script_path', record)
        self.assertIn('retrieved_utc', record)
        self.assertNotIn('metadata_kind', record)
        self.assertEqual(json.loads(self.mapping.provenance_path(output).read_text()), record)
        self.assertEqual(list(self.work.glob('*.acquisition.*.json')), [])

    def test_existing_canonical_service_legacy_receipt_is_not_sanitized_or_rewritten(self):
        output, raw, _, constants = self.service('hg19')
        full, = self.work.glob(output.name + '.acquisition.*.json')
        metadata = self.mapping.provenance_path(output)
        metadata.write_bytes(full.read_bytes())
        before = metadata.read_bytes(), metadata.stat().st_mtime_ns, output.read_bytes()
        unused_legacy_raw = output.with_name('historical-unused.raw.tsv')
        with mock.patch.multiple(self.mapping, **constants), \
                mock.patch.object(self.mapping, '_request', side_effect=AssertionError('no reacquisition')):
            self.mapping.acquire(output, unused_legacy_raw)
        self.assertEqual(before, (metadata.read_bytes(), metadata.stat().st_mtime_ns, output.read_bytes()))
        self.assertEqual(len(list(self.work.glob('*.acquisition.*.json'))), 1)
        self.assertFalse(unused_legacy_raw.exists())

    def test_new_canonical_ftp_projection_has_portable_metadata_and_external_full_receipt(self):
        output, record, _ = self.ftp()
        self.assert_portable(record)
        self.assertIn('source_root', record)
        self.assertIn('raw_sources', record)
        self.assertIn('schema_sha256', record)
        self.assertIn('query_sql', record)
        self.assertIs(record['online_query_executed'], False)
        self.assertIs(record['original_biomart_response_order_identity'], False)
        full, = self.work.glob(output.name + '.acquisition.*.json')
        runtime = json.loads(full.read_text())
        self.assertIn('registration_producer_path', runtime)
        self.assertIn('registered_utc', runtime)
        self.assertIn('slurm_job_id', runtime)
        self.assertEqual(output.read_bytes(), BODY)

    def test_new_external_ftp_projection_keeps_its_original_full_sidecar_behavior(self):
        output, record, _ = self.ftp(canonical=False)
        self.assertIn('registration_producer_path', record)
        self.assertIn('registered_utc', record)
        self.assertNotIn('metadata_kind', record)
        self.assertEqual(json.loads(self.mapping.provenance_path(output).read_text()), record)
        self.assertEqual(list(self.work.glob('*.acquisition.*.json')), [])

    def test_existing_canonical_ftp_legacy_receipt_is_not_sanitized_or_rewritten(self):
        output, _, constants = self.ftp()
        full, = self.work.glob(output.name + '.acquisition.*.json')
        metadata = self.mapping.provenance_path(output)
        metadata.write_bytes(full.read_bytes())
        before = metadata.read_bytes(), metadata.stat().st_mtime_ns, output.read_bytes()
        unused_legacy_workspace = self.checkout / 'historical-unused-workspace'
        with mock.patch.multiple(self.mapping.FTP, **constants), \
                mock.patch.object(self.mapping.FTP, 'acquire', side_effect=AssertionError('no reacquisition')):
            self.mapping.acquire(output, source_id=self.mapping.FTP.SOURCE_ID, resource_dir=unused_legacy_workspace)
        self.assertEqual(before, (metadata.read_bytes(), metadata.stat().st_mtime_ns, output.read_bytes()))
        self.assertEqual(len(list(self.work.glob('*.acquisition.*.json'))), 1)
        self.assertFalse(unused_legacy_workspace.exists())

    def test_external_runtime_receipts_never_replace_an_earlier_same_basename_record(self):
        first = self.mapping.runtime_receipt({'run': 'first'}, self.work, 'same.tsv')
        before = first.read_bytes(), first.stat().st_mtime_ns
        second = self.mapping.runtime_receipt({'run': 'second'}, self.work, 'same.tsv')
        self.assertNotEqual(first, second)
        self.assertEqual(before, (first.read_bytes(), first.stat().st_mtime_ns))
        self.assertEqual(json.loads(second.read_text()), {'run': 'second'})


if __name__ == '__main__':
    unittest.main()
