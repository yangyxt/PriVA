"""Exercise ClinVar installer stage order and failure publication with real headers/pickles."""
from __future__ import annotations

import gzip
import hashlib
import json
import os
from pathlib import Path
import pickle
import subprocess
import sys
import tempfile
import unittest
import yaml

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'tests'))
import fake_priva_binary
from test_clinvar_regeneration import header_vcf, row
from clinvar_annotation_contract import AA_CSQ_FIELDS, GENE_CSQ_FIELDS, AF_CSQ_FIELDS
from clinvar_stat_variants import per_gene_stat, logger
import pandas as pd
import pysam


PRELUDE = r'''
set -euo pipefail
event() { printf '%s\n' "$*" >> "$TEST_EVENTS"; }
log() { printf '%s\n' "$*"; }
randomID() { printf 'fixture'; }
read_yaml() {
    case "$2" in
        assembly) printf hg19 ;;
        clinvar_vcf_dir) printf '%s' "$TEST_CACHE" ;;
        clinvar_vcf) printf '%s' "$TEST_VCF" ;;
        clinvar_aa_stat) printf '%s' "$TEST_AA" ;;
        clinvar_splice_stat) printf '%s' "$TEST_SPLICE" ;;
        clinvar_gene_stat) printf '%s' "$TEST_BASELINE_CONFIG" ;;
        clinvar_patho_af_stat) printf '%s' "${TEST_GENE_AF:-null}" ;;
        clinvar_patho_exon_af_stat) printf '%s' "${TEST_EXON_AF:-null}" ;;
        priva_python_hash_seed) printf '%s' "${TEST_SEED:-null}" ;;
        threads) printf 1 ;;
        gnomad_vcf_chrX) printf '%s/gnomad.chrX.vcf.gz' "$TEST_CACHE" ;;
        *) printf null ;;
    esac
}
update_yaml() { event "update $2 $3"; }
# The fixture priva (tests/fake_priva_binary.py) records each native call and
# runs the replaced Python script through the exported python function below.
priva_resolve_binary() { PRIVA_BIN=$PRIVA_BINARY; }
check_vcf_validity() { [[ -s "$1" && -s "$1.tbi" ]]; }
download_resource_atomic() {
    event "download $1"
    if [[ "$2" == *.tbi ]]; then cp "$TEST_RAW.tbi" "$2"; else cp "$TEST_RAW" "$2"; fi
}
liftover_from_ucsc_to_GRCh() { event rename; cp "$1" "$3"; cp "$1.tbi" "$3.tbi"; }
bash() {
    case "$2" in
        anno_VEP_cache_candidate) event vep; cp "$TEST_VEP" "$4"; cp "$TEST_VEP.tbi" "$4.tbi" ;;
        anno_agg_gnomAD_data)
            event gnomad
            [[ "$TEST_FAIL" != gnomad ]] || return 7
            cp "$TEST_COMPLETE" "$3"; cp "$TEST_COMPLETE.tbi" "$3.tbi" ;;
        *) return 8 ;;
    esac
}
python() {
    case "${1##*/}" in
        clinvar_annotation_contract.py) "$TEST_PYTHON" "$@" ;;
        clinvar_gene_refresh.py) "$TEST_PYTHON" "$@" ;;
        stat_aachange_clinvar.py)
            event aa_builder
            [[ "$TEST_FAIL" != aa ]] || return 7
            cp "$TEST_AA_FIXTURE" "$3"; cp "$TEST_SPLICE_FIXTURE" "$4" ;;
        stat_gene_patho_afs.py)
            event gene_af_builder
            # Without a CSQ format the builder succeeds and writes nothing.
            [[ "$TEST_FAIL" != gene_af_no_csq ]] || return 0
            cp "$TEST_GENE_AF_FIXTURE" "$3" ;;
        stat_exon_patho_afs.py)
            event exon_af_builder
            [[ "$TEST_FAIL" != exon_af ]] || return 7
            cp "$TEST_EXON_AF_FIXTURE" "$5" ;;
        clinvar_stat_variants.py)
            event "gene_builder $*"
            local out= argument
            while [[ $# -gt 0 ]]; do
                argument=$1; shift
                if [[ "$argument" == --output_pickle ]]; then out=$1; fi
            done
            cp "$TEST_GENE_NATIVE_FIXTURE" "$out" ;;
        *) return 9 ;;
    esac
}
export -f event python
'''


class ClinVarInstallerTests(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory()
        self.addCleanup(self.temporary.cleanup)
        self.directory = Path(self.temporary.name)
        self.cache = self.directory / 'cache'
        self.cache.mkdir()
        self.events = self.directory / 'events.txt'
        self.events.touch()
        fields = AA_CSQ_FIELDS | GENE_CSQ_FIELDS | AF_CSQ_FIELDS
        paths = {}
        for name, csq, population in [('raw', {'Feature', 'SYMBOL'}, False),
                                      ('vep', fields, False), ('complete', fields, True)]:
            plain = self.directory / (name + '.vcf')
            header_vcf(plain, csq, population=population)
            if name != 'raw':
                plain.write_text(plain.read_text().replace('##fileformat=VCFv4.2\n',
                    '##fileformat=VCFv4.2\n##VEP="v116.3" assembly="GRCh37"\n'))
            path = self.directory / (name + '.vcf.gz')
            pysam.tabix_compress(str(plain), str(path))
            pysam.tabix_index(str(path), preset='vcf')
            paths[name] = path
        data = self.directory / 'data'
        (data / 'ClinVar').mkdir(parents=True)
        (data / 'liftover').mkdir()
        (data / 'liftover/GRC_to_ucsc.contig.map.txt').write_text('1\tchr1\n')
        self.baseline = data / 'ClinVar/clinvar_2star_stats.pkl.gz'
        logger.disabled = True
        gene = per_gene_stat(pd.DataFrame([row()]))
        fixtures = {
            'aa': {'ENST1': {'1': {'p.A1V': {'CLNSIG': ['Pathogenic'], 'CLNREVSTAT': ['reviewed_by_expert_panel']}}}},
            'splice': {'ENST1': [{'chrom': 'chr1', 'pos': 1, 'clinvar_sig': 'Pathogenic', 'splice_ai': {}}]},
            'gene_af': {'ENSG1': None},
            'exon_af': {'ENST1': {'1/2': (0., 0., 0., 0., [])}},
            'gene': {'ENSG00000000001': gene},
        }
        fixture_paths = {}
        for name, value in fixtures.items():
            path = self.directory / (name + ('.pkl.gz' if name == 'gene' else '.pkl'))
            opener = gzip.open if name == 'gene' else open
            with opener(path, 'wb') as handle:
                pickle.dump(value, handle)
            fixture_paths[name] = path
        # The fixture priva checks its archives with the Python cache checks,
        # which read a plain pickle at a .rkyv name.
        fixture_paths['gene_native'] = self.directory / 'gene_native.pkl'
        fixture_paths['gene_native'].write_bytes(pickle.dumps(fixtures['gene']))
        self.baseline.write_bytes(fixture_paths['gene'].read_bytes())
        priva = tempfile.TemporaryDirectory()
        self.addCleanup(priva.cleanup)
        self.env = dict(os.environ, TEST_EVENTS=str(self.events), TEST_CACHE=str(self.cache),
                        TEST_VCF=str(self.cache / 'clinvar.hg19.vep.vcf.gz'),
                        TEST_AA=str(self.cache / 'clinvar.hg19.vep.aa_change.stats.rkyv'),
                        TEST_SPLICE=str(self.cache / 'clinvar.hg19.vep.splice_change.stats.rkyv'),
                        TEST_BASELINE_CONFIG='clinvar_2star_stats.pkl.gz',
                        TEST_RAW=str(paths['raw']), TEST_VEP=str(paths['vep']),
                        TEST_COMPLETE=str(paths['complete']), TEST_PYTHON=sys.executable,
                        TEST_FAIL='none', SCRIPT_DIR=str(ROOT / 'scripts'), DATA_DIR=str(data), BASE_DIR=str(ROOT),
                        **fake_priva_binary.environment(priva.name))
        self.env.update({'TEST_' + name.upper() + '_FIXTURE': str(path)
                         for name, path in fixture_paths.items()})
        source = (ROOT / 'scripts/install_utils.sh').read_text()
        self.functions = source[source.index('function ClinVar_VCF_deploy() {'):
                                source.index('function Conservation_install()')]
        self.functions += source[source.index('function yaml_value_or_default() {'):
                                 source.index('function update_or_append_yaml() {')]

    def run_function(self, command, fail='none'):
        # These existing fixtures isolate the scientific body and its stage
        # order. Public ownership/retirement routing has its own real fixtures.
        if command.startswith('ClinVar_VCF_deploy '):
            config = self.directory / 'deployment.yaml'
            config.write_text(yaml.safe_dump(dict(cache_root=str(self.cache),
                clinvar_vcf_dir=str(self.cache), clinvar_vcf=self.env['TEST_VCF'],
                assembly='hg19', vep_installed_version='116')))
            command = command.replace('ClinVar_VCF_deploy config.yaml',
                                      '_ClinVar_VCF_deploy "' + str(config) + '"', 1)
        environment = dict(self.env, TEST_FAIL=fail)
        result = subprocess.run(['bash', '-c', PRELUDE + self.functions + '\n' + command],
                                env=environment, text=True, capture_output=True)
        return result, self.events.read_text().splitlines()

    def priva_calls(self):
        calls = fake_priva_binary.calls(self.env['PRIVA_FAKE_CALLS'])
        Path(self.env['PRIVA_FAKE_CALLS']).unlink(missing_ok=True)
        return calls

    def contract_calls(self):
        """priva calls, with each AA/splice pair call cut to its subcommand.

        Pair calls must name the Python contract, whose SHA-256 the receipts keep.
        """
        calls = self.priva_calls()
        for call in calls:
            if call[:1] == ['clinvar-contract'] and call[1].startswith('pair-'):
                source = call[call.index('--contract-source') + 1]
                self.assertEqual(Path(source).name, 'clinvar_annotation_contract.py')
        return [call[:2] if call[:1] == ['clinvar-contract'] and call[1].startswith('pair-') else call
                for call in calls]

    def install_vcf(self, stage):
        path = Path(self.env['TEST_' + stage.upper()])
        target = Path(self.env['TEST_VCF'])
        target.write_bytes(path.read_bytes())
        Path(str(target) + '.tbi').write_bytes(Path(str(path) + '.tbi').read_bytes())

    def test_fresh_vcf_preserves_download_rename_full_vep_gnomad_order(self):
        result, events = self.run_function('ClinVar_VCF_deploy config.yaml')
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertEqual([event.split()[0] for event in events],
                         ['download', 'download', 'rename', 'vep', 'gnomad', 'update'])

    def test_pinned_source_uses_original_annotation_order_without_moving_download_url(self):
        result, events = self.run_function('ClinVar_VCF_deploy config.yaml --source-vcf "$TEST_RAW"')
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertEqual([event.split()[0] for event in events], ['rename', 'vep', 'gnomad', 'update'])

    def test_pinned_source_cannot_replace_an_explicit_annotated_vcf(self):
        custom = self.directory / 'selected.vcf.gz'
        custom.write_bytes(Path(self.env['TEST_COMPLETE']).read_bytes())
        self.env['TEST_VCF'] = str(custom)
        before = custom.read_bytes()
        result, events = self.run_function('ClinVar_VCF_deploy config.yaml --source-vcf "$TEST_RAW"')
        self.assertNotEqual(result.returncode, 0)
        self.assertIn('cannot replace an explicitly selected', result.stdout)
        self.assertEqual(custom.read_bytes(), before)
        self.assertEqual(events, [])

    def test_complete_vep_stage_resumes_at_gnomad(self):
        self.install_vcf('vep')
        result, events = self.run_function('ClinVar_VCF_deploy config.yaml')
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertEqual([event.split()[0] for event in events], ['gnomad', 'update'])

    def test_complete_default_cache_registers_resolved_path(self):
        self.install_vcf('complete')
        actual = self.env['TEST_VCF']
        self.env['TEST_VCF'] = Path(actual).name
        result, events = self.run_function('ClinVar_VCF_deploy config.yaml')
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertEqual(events, ['update clinvar_vcf ' + actual])

    def test_explicit_complete_input_is_honored_without_cache_substitution(self):
        self.install_vcf('complete')
        custom = self.directory / Path(self.env['TEST_VCF']).name
        custom.write_bytes(Path(self.env['TEST_COMPLETE']).read_bytes())
        self.env['TEST_VCF'] = str(custom)
        before = custom.read_bytes()
        result, events = self.run_function('ClinVar_VCF_deploy config.yaml')
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertEqual(events, [])
        self.assertEqual(custom.read_bytes(), before)
        self.assertIn(str(custom), result.stdout)

    def test_explicit_missing_or_incomplete_input_cannot_fall_back_to_complete_cache(self):
        self.install_vcf('complete')
        custom = self.directory / 'custom.vcf.gz'
        self.env['TEST_VCF'] = str(custom)
        for present in (False, True):
            with self.subTest(present=present):
                if present:
                    custom.write_bytes(Path(self.env['TEST_RAW']).read_bytes())
                before = custom.read_bytes() if custom.exists() else None
                result, events = self.run_function('ClinVar_VCF_deploy config.yaml')
                self.assertNotEqual(result.returncode, 0)
                self.assertEqual(events, [])
                self.assertEqual(custom.read_bytes() if custom.exists() else None, before)

    def test_noncanonical_basename_is_explicit_and_not_overwritten(self):
        custom = self.cache / 'custom.vcf.gz'
        custom.write_bytes(Path(self.env['TEST_RAW']).read_bytes())
        self.env['TEST_VCF'] = custom.name
        before = custom.read_bytes()
        result, events = self.run_function('ClinVar_VCF_deploy config.yaml')
        self.assertNotEqual(result.returncode, 0)
        self.assertEqual(events, [])
        self.assertEqual(custom.read_bytes(), before)

    def test_explicit_corrupt_gzip_cannot_be_repaired_from_plain_sibling(self):
        custom = self.directory / 'custom.vcf.gz'
        custom.write_bytes(b'corrupt gzip retained')
        Path(str(custom)[:-3]).write_text('newer plain sibling must not replace custom bytes\n')
        self.env['TEST_VCF'] = str(custom)
        result, events = self.run_function('ClinVar_VCF_deploy config.yaml')
        self.assertNotEqual(result.returncode, 0)
        self.assertEqual(events, [])
        self.assertEqual(custom.read_bytes(), b'corrupt gzip retained')

    def test_gnomad_failure_blocks_configuration_success(self):
        self.install_vcf('vep')
        result, events = self.run_function('ClinVar_VCF_deploy config.yaml', fail='gnomad')
        self.assertNotEqual(result.returncode, 0)
        self.assertEqual(events, ['gnomad'])

    def test_aa_writer_failure_preserves_previous_outputs(self):
        self.install_vcf('complete')
        for key in ('TEST_AA', 'TEST_SPLICE'):
            Path(self.env[key]).write_bytes(b'previous')
        result, events = self.run_function('ClinVar_AA_stat config.yaml', fail='aa')
        self.assertNotEqual(result.returncode, 0)
        self.assertEqual(events, ['aa_builder'])
        self.assertEqual(Path(self.env['TEST_AA']).read_bytes(), b'previous')
        self.assertEqual(Path(self.env['TEST_SPLICE']).read_bytes(), b'previous')

    def test_valid_aa_and_splice_publish_at_established_names(self):
        self.install_vcf('complete')
        result, events = self.run_function('ClinVar_AA_stat config.yaml')
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertEqual([event.split()[0] for event in events], ['aa_builder', 'update', 'update'])
        self.assertTrue(Path(self.env['TEST_AA']).is_file())
        self.assertTrue(Path(self.env['TEST_SPLICE']).is_file())
        vcf, aa, splice = self.env['TEST_VCF'], self.env['TEST_AA'], self.env['TEST_SPLICE']
        self.assertEqual(self.contract_calls(), [
            ['clinvar-contract', 'header', '--vcf', vcf, '--stage', 'aa'],
            ['builder-identity', 'clinvar-aa-splice'],
            ['clinvar-contract', 'pair-capture'],
            ['clinvar-aa-splice', 'build', '--vcf', vcf, '--aa-output', aa + '.partial.fixture',
             '--splice-output', splice + '.partial.fixture', '--threads', '1'],
            ['clinvar-contract', 'cache', '--input', aa + '.partial.fixture', '--kind', 'aa'],
            ['clinvar-contract', 'cache', '--input', splice + '.partial.fixture', '--kind', 'splice'],
            ['clinvar-contract', 'pair-complete']])
        receipt = json.loads(Path(aa + '.pair.ready.json').read_text())
        producer = Path(aa + '.builder_identity.json')
        self.assertEqual(receipt['producer_path'], str(producer))
        self.assertEqual(json.loads(producer.read_text())['builder'], 'clinvar-aa-splice')
        self.assertEqual(receipt['inputs']['producer_sha256'], hashlib.sha256(producer.read_bytes()).hexdigest())

    def test_native_aa_reuse_follows_the_builder_identity(self):
        self.install_vcf('complete')
        result, _ = self.run_function('ClinVar_AA_stat config.yaml')
        self.assertEqual(result.returncode, 0, result.stderr)
        self.priva_calls()
        self.events.write_text('')
        result, events = self.run_function('ClinVar_AA_stat config.yaml')
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertEqual(events, [])
        vcf, aa, splice = self.env['TEST_VCF'], self.env['TEST_AA'], self.env['TEST_SPLICE']
        self.assertEqual(self.contract_calls(), [
            ['clinvar-contract', 'header', '--vcf', vcf, '--stage', 'aa'],
            ['builder-identity', 'clinvar-aa-splice'],
            ['clinvar-contract', 'cache', '--input', aa, '--kind', 'aa'],
            ['clinvar-contract', 'cache', '--input', splice, '--kind', 'splice'],
            ['clinvar-contract', 'pair-check']])
        self.env['PRIVA_FAKE_BUILDER_SHA256'] = 'f' * 64
        result, events = self.run_function('ClinVar_AA_stat config.yaml')
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertIn('aa_builder', events)

    def test_legacy_python_pair_is_reused_without_native_builds(self):
        self.install_vcf('complete')
        aa, splice = self.cache / 'legacy.aa.pkl', self.cache / 'legacy.splice.pkl'
        for path, fixture in ((aa, 'TEST_AA_FIXTURE'), (splice, 'TEST_SPLICE_FIXTURE')):
            path.write_bytes(Path(self.env[fixture]).read_bytes())
        contract = [sys.executable, str(ROOT / 'scripts/clinvar_annotation_contract.py')]
        pair = ['--vcf', self.env['TEST_VCF'], '--producer', str(ROOT / 'scripts/stat_aachange_clinvar.py'),
                '--aa', str(aa), '--splice', str(splice), '--receipt', str(aa) + '.pair.ready.json',
                '--captured', str(self.directory / 'captured.json')]
        for operation in ('pair-capture', 'pair-complete'):
            subprocess.run(contract + [operation] + pair, check=True, capture_output=True)
        self.env.update(TEST_AA=str(aa), TEST_SPLICE=str(splice))
        result, events = self.run_function('ClinVar_AA_stat config.yaml')
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertEqual(events, [])
        self.assertEqual(self.contract_calls(), [
            ['clinvar-contract', 'header', '--vcf', self.env['TEST_VCF'], '--stage', 'aa'],
            ['clinvar-contract', 'legacy-cache', '--input', str(aa), '--kind', 'aa'],
            ['clinvar-contract', 'legacy-cache', '--input', str(splice), '--kind', 'splice'],
            ['clinvar-contract', 'pair-check']])
        self.assertFalse(Path(str(aa)[:-4] + '.rkyv').exists())

    def test_aa_outputs_honor_distinct_configured_directories(self):
        self.install_vcf('complete')
        # Pickle values restored by cache-install prepare select .rkyv siblings.
        self.env['TEST_AA'] = str(self.cache / 'aa/custom.pkl')
        self.env['TEST_SPLICE'] = str(self.cache / 'splice/custom.pkl')
        result, events = self.run_function('ClinVar_AA_stat config.yaml')
        self.assertEqual(result.returncode, 0, result.stderr)
        aa, splice = self.cache / 'aa/custom.rkyv', self.cache / 'splice/custom.rkyv'
        self.assertEqual(events, ['aa_builder', 'update clinvar_aa_stat ' + str(aa),
                                  'update clinvar_splice_stat ' + str(splice)])
        self.assertEqual(aa.read_bytes(), Path(self.env['TEST_AA_FIXTURE']).read_bytes())
        self.assertEqual(splice.read_bytes(), Path(self.env['TEST_SPLICE_FIXTURE']).read_bytes())
        self.assertFalse(Path(self.env['TEST_AA']).exists())
        self.assertFalse((self.cache / 'clinvar.hg19.vep.aa_change.stats.rkyv').exists())

    def test_managed_aa_archive_is_never_rebuilt_in_place(self):
        self.install_vcf('complete')
        managed = {}
        for key in ('TEST_AA', 'TEST_SPLICE'):
            managed[key] = self.cache / ('managed.' + key.lower() + '.v1.rkyv')
            managed[key].write_bytes(b'published archive')
            Path(str(managed[key]) + '.source.json').write_text('{"origin": "direct"}\n')
            self.env[key] = str(managed[key])
        result, events = self.run_function('ClinVar_AA_stat config.yaml')
        self.assertEqual(result.returncode, 0, result.stderr)
        defaults = [self.cache / 'clinvar.hg19.vep.aa_change.stats.rkyv',
                    self.cache / 'clinvar.hg19.vep.splice_change.stats.rkyv']
        self.assertEqual(events, ['aa_builder', 'update clinvar_aa_stat ' + str(defaults[0]),
                                  'update clinvar_splice_stat ' + str(defaults[1])])
        self.assertTrue(all(path.is_file() for path in defaults))
        self.assertTrue(all(path.read_bytes() == b'published archive' for path in managed.values()))

    def test_aa_pair_completion_cold_warm_and_input_output_invalidation(self):
        self.install_vcf('complete')
        for iteration, changed in enumerate((None, None, 'TEST_VCF', 'TEST_AA', 'TEST_SPLICE')):
            if changed:
                path = Path(self.env[changed])
                os.utime(path, ns=(path.stat().st_atime_ns, path.stat().st_mtime_ns + 1000000))
            self.events.write_text('')
            result, events = self.run_function('ClinVar_AA_stat config.yaml')
            self.assertEqual(result.returncode, 0, result.stderr)
            self.assertEqual('aa_builder' in events, iteration != 1)
            for key in ('AA', 'SPLICE'):
                self.assertEqual(Path(self.env['TEST_' + key]).read_bytes(),
                                 Path(self.env['TEST_' + key + '_FIXTURE']).read_bytes())

    def test_aa_pair_old_schema_valid_outputs_require_latest_builder(self):
        self.install_vcf('complete')
        for key in ('AA', 'SPLICE'):
            Path(self.env['TEST_' + key]).write_bytes(Path(self.env['TEST_' + key + '_FIXTURE']).read_bytes())
        result, events = self.run_function('ClinVar_AA_stat config.yaml', fail='aa')
        self.assertNotEqual(result.returncode, 0)
        self.assertIn('aa_builder', events)
        self.assertFalse(Path(self.env['TEST_AA'] + '.pair.ready.json').exists())

    def test_aa_pair_yaml_failure_has_no_completion_then_retry(self):
        self.install_vcf('complete')
        result, _ = self.run_function('update_yaml() { [[ "$2" != clinvar_splice_stat ]]; }; ClinVar_AA_stat config.yaml')
        self.assertNotEqual(result.returncode, 0)
        self.assertFalse(Path(self.env['TEST_AA'] + '.pair.ready.json').exists())
        result, events = self.run_function('ClinVar_AA_stat config.yaml')
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertIn('aa_builder', events)

    def test_aa_pair_content_identity_cross_checkout_and_changed_producer(self):
        from clinvar_annotation_contract import pair_inputs
        self.install_vcf('complete')
        producer = ROOT / 'scripts/stat_aachange_clinvar.py'
        copied = self.directory / 'copied.py'
        copied.write_bytes(producer.read_bytes())
        original = pair_inputs(self.env['TEST_VCF'], str(producer))
        self.assertEqual(original, pair_inputs(self.env['TEST_VCF'], str(copied)))
        copied.write_bytes(copied.read_bytes() + b'\n# changed producer\n')
        self.assertNotEqual(original, pair_inputs(self.env['TEST_VCF'], str(copied)))
        index = Path(self.env['TEST_VCF'] + '.tbi')
        os.utime(index, ns=(index.stat().st_atime_ns, index.stat().st_mtime_ns + 1))
        self.assertNotEqual(original, pair_inputs(self.env['TEST_VCF'], str(producer)))

    def test_aa_output_collision_stops_before_builder(self):
        self.install_vcf('complete')
        self.env['TEST_AA'] = self.env['TEST_SPLICE']
        result, events = self.run_function('ClinVar_AA_stat config.yaml')
        self.assertNotEqual(result.returncode, 0)
        self.assertEqual(events, [])

    def test_aa_receipt_named_splice_value_writes_its_own_archive(self):
        self.install_vcf('complete')
        receipt = Path(self.env['TEST_AA'] + '.pair.ready.json')
        self.env['TEST_SPLICE'] = str(receipt)
        result, events = self.run_function('ClinVar_AA_stat config.yaml')
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertIn('aa_builder', events)
        self.assertEqual(json.loads(receipt.read_text())['schema'], 1)
        self.assertEqual(Path(str(receipt) + '.rkyv').read_bytes(),
                         Path(self.env['TEST_SPLICE_FIXTURE']).read_bytes())

    def test_aa_linked_producer_file_is_replaced_not_written_through(self):
        self.install_vcf('complete')
        producer = Path(self.env['TEST_AA'] + '.builder_identity.json')
        protected = self.directory / 'protected.txt'
        protected.write_text('linked file remains unchanged\n')
        os.symlink(protected, producer)
        result, events = self.run_function('ClinVar_AA_stat config.yaml')
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertIn('aa_builder', events)
        self.assertFalse(producer.is_symlink())
        self.assertEqual(json.loads(producer.read_text())['builder'], 'clinvar-aa-splice')
        self.assertEqual(protected.read_text(), 'linked file remains unchanged\n')

    def test_aa_real_producer_pair_bytes_order_and_warm_reuse(self):
        fields = sorted(AA_CSQ_FIELDS) + ['Unused']
        values = dict(Feature='ENST1', Feature_type='Transcript', Protein_position='10',
                      HGVSp='p.A10C', Consequence='missense_variant', EXON='1/2', HGVSc='c.1A>C')
        header = ('##fileformat=VCFv4.2\n##contig=<ID=1>\n'
                  '##INFO=<ID=CLNSIG,Number=.,Type=String,Description="fixture">\n'
                  '##INFO=<ID=CLNREVSTAT,Number=.,Type=String,Description="fixture">\n'
                  '##INFO=<ID=CSQ,Number=.,Type=String,Description="Format: ' + '|'.join(fields) + '">\n'
                  '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n')
        raw = self.directory / 'real.vcf'
        raw.write_text(header + ''.join('1\t' + str(pos) + '\t.\tA\tC\t.\tPASS\tCLNSIG=Pathogenic;'
                       'CLNREVSTAT=reviewed_by_expert_panel;CSQ=' + '|'.join(values.get(f, '') for f in fields) + '\n'
                       for pos in (10, 20)))
        pysam.tabix_compress(str(raw), self.env['TEST_VCF'], force=True)
        pysam.tabix_index(self.env['TEST_VCF'], preset='vcf', force=True)
        expected = [self.directory / 'expected-aa.pkl', self.directory / 'expected-splice.pkl']
        baseline = Path(os.environ.get('PRIVA_WORKER_BASELINE', str(ROOT)))
        subprocess.run([sys.executable, str(baseline / 'scripts/stat_aachange_clinvar.py'),
                        self.env['TEST_VCF'], *map(str, expected), '1'], check=True, capture_output=True)
        command = '''python() { [[ "${1##*/}" != stat_aachange_clinvar.py ]] || event aa_builder; "$TEST_PYTHON" "$@"; }; export -f python; ClinVar_AA_stat config.yaml'''
        result, events = self.run_function(command)
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertIn('aa_builder', events)
        for key, previous in zip(('TEST_AA', 'TEST_SPLICE'), expected):
            self.assertEqual(Path(self.env[key]).read_bytes(), previous.read_bytes())
        self.events.write_text('')
        result, events = self.run_function(command)
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertNotIn('aa_builder', events)

    def test_aa_pair_config_and_hardlink_aliases_fail_without_mutation(self):
        self.install_vcf('complete')
        receipt = Path(self.env['TEST_AA'] + '.pair.ready.json')
        config = self.directory / 'selected.yaml'
        config.write_text('original YAML\n')
        for kind in ('same-path', 'hardlink', 'symlink'):
            with self.subTest(kind=kind):
                if receipt.exists() or receipt.is_symlink():
                    receipt.unlink()
                if kind == 'same-path':
                    receipt.write_bytes(config.read_bytes())
                    selected = receipt
                else:
                    (os.link if kind == 'hardlink' else os.symlink)(config, receipt)
                    selected = config
                self.events.write_text('')
                result, events = self.run_function('ClinVar_AA_stat "' + str(selected) + '"')
                self.assertNotEqual(result.returncode, 0)
                self.assertEqual(events, [])
                self.assertEqual(config.read_text(), 'original YAML\n')
                self.assertEqual(receipt.read_text(), 'original YAML\n')

    def test_aa_second_move_failure_has_no_completion_and_retry_rebuilds(self):
        self.install_vcf('complete')
        command = '''mv() { [[ "${@: -1}" != "$TEST_SPLICE" ]] || return 7; command mv "$@"; }; ClinVar_AA_stat config.yaml'''
        result, _ = self.run_function(command)
        self.assertNotEqual(result.returncode, 0)
        self.assertFalse(Path(self.env['TEST_AA'] + '.pair.ready.json').exists())
        result, events = self.run_function('ClinVar_AA_stat config.yaml')
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertIn('aa_builder', events)

    def af_case(self):
        self.install_vcf('complete')
        config = self.directory / 'active config.yaml'
        config.write_text('active configuration remains unchanged\n')
        self.env['TEST_CONFIG'] = str(config)
        return config

    def fixture_state(self):
        result = {}
        for path in [self.directory, *sorted(self.directory.rglob('*'))]:
            info = path.lstat()
            physical = tuple(getattr(info, 'st_' + key) for key in
                             ('dev', 'ino', 'mode', 'nlink', 'size', 'mtime_ns', 'ctime_ns'))
            content = os.readlink(path) if path.is_symlink() else path.read_bytes() if path.is_file() else None
            result[str(path.relative_to(self.directory))] = (physical, content)
        return result

    def assert_af_alias_rejected(self, output, target_name, kind):
        # The builders write .rkyv archives, so the aliases are at those names.
        case = self.cache / ('af-case-' + output + '-' + target_name + '-' + kind)
        gene, exon = case / 'gene/result.rkyv', case / 'exon/result.rkyv'
        other = exon if output == 'gene' else gene
        target = {'vcf': Path(self.env['TEST_VCF']),
                  'tbi': Path(self.env['TEST_VCF'] + '.tbi'),
                  'csi': Path(self.env['TEST_VCF'] + '.csi'),
                  'config': Path(self.env['TEST_CONFIG']), 'other': other}[target_name]
        if kind in ('symlink', 'hardlink') and not target.exists():
            target.parent.mkdir(parents=True, exist_ok=True)
            target.write_bytes(b'pre-existing protected file\n')
        selected = gene if output == 'gene' else exon
        if kind == 'normalized':
            selected = target.parent / 'must-not-be-created' / '..' / target.name
        elif kind == 'same':
            selected = target
        else:
            selected.parent.mkdir(parents=True, exist_ok=True)
            (os.symlink if kind == 'symlink' else os.link)(target, selected)
        self.env.update(TEST_GENE_AF=str(selected if output == 'gene' else gene),
                        TEST_EXON_AF=str(selected if output == 'exon' else exon))
        before = self.fixture_state()
        result, events = self.run_function('ClinVar_patho_AF_stat "$TEST_CONFIG"')
        self.assertNotEqual(result.returncode, 0)
        self.assertRegex(result.stderr, 'paths must (be distinct|not share an existing file)')
        self.assertEqual(events, [])
        self.assertEqual(self.fixture_state(), before)
        self.assertEqual(self.priva_calls(), [
            ['clinvar-contract', 'af-paths-check', '--vcf', self.env['TEST_VCF'], '--config', self.env['TEST_CONFIG'],
             '--gene', self.env['TEST_GENE_AF'], '--exon', self.env['TEST_EXON_AF']]])

    def test_af_normalized_output_aliases_fail_without_writes(self):
        self.af_case()
        for output in ('gene', 'exon'):
            with self.subTest(output=output):
                self.assert_af_alias_rejected(output, 'other', 'normalized')

    def test_af_symlink_input_index_config_and_output_aliases_fail_without_writes(self):
        self.af_case()
        for output in ('gene', 'exon'):
            for target in ('vcf', 'tbi', 'csi', 'config', 'other'):
                with self.subTest(output=output, target=target):
                    self.assert_af_alias_rejected(output, target, 'symlink')

    def test_af_hardlink_input_index_config_and_output_aliases_fail_without_writes(self):
        self.af_case()
        for output in ('gene', 'exon'):
            for target in ('vcf', 'tbi', 'csi', 'config', 'other'):
                with self.subTest(output=output, target=target):
                    self.assert_af_alias_rejected(output, target, 'hardlink')

    def test_af_same_missing_output_path_is_rejected_before_directory_creation(self):
        self.af_case()
        for output in ('gene', 'exon'):
            with self.subTest(output=output):
                self.assert_af_alias_rejected(output, 'other', 'same')

    def test_af_values_naming_protected_files_write_separate_archives(self):
        config = self.af_case()
        protected = [config, Path(self.env['TEST_VCF']), Path(self.env['TEST_VCF'] + '.tbi')]
        before = {path: path.read_bytes() for path in protected}
        self.env.update(TEST_GENE_AF=self.env['TEST_VCF'], TEST_EXON_AF=str(config))
        result, events = self.run_function('ClinVar_patho_AF_stat "$TEST_CONFIG"')
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertEqual(events, ['gene_af_builder', 'exon_af_builder',
                                  'update clinvar_patho_af_stat ' + self.env['TEST_VCF'] + '.rkyv',
                                  'update clinvar_patho_exon_af_stat ' + str(config) + '.rkyv'])
        self.assertEqual({path: path.read_bytes() for path in protected}, before)

    def test_af_distinct_configured_paths_keep_original_bytes_and_producer_order(self):
        config = self.af_case()
        protected = [config, Path(self.env['TEST_VCF']), Path(self.env['TEST_VCF'] + '.tbi')]
        before = {path: path.read_bytes() for path in protected}
        # Pickle values restored by cache-install prepare select .rkyv siblings.
        gene, exon = self.cache / "gene's output/result.pkl", self.cache / 'exon output/result.pkl'
        self.env.update(TEST_GENE_AF=str(gene), TEST_EXON_AF=str(exon))
        result, events = self.run_function('ClinVar_patho_AF_stat "$TEST_CONFIG"')
        self.assertEqual(result.returncode, 0, result.stderr)
        gene, exon = gene.with_suffix('.rkyv'), exon.with_suffix('.rkyv')
        self.assertEqual(events, ['gene_af_builder', 'exon_af_builder',
                                  'update clinvar_patho_af_stat ' + str(gene),
                                  'update clinvar_patho_exon_af_stat ' + str(exon)])
        for output, fixture in ((gene, 'TEST_GENE_AF_FIXTURE'), (exon, 'TEST_EXON_AF_FIXTURE')):
            self.assertEqual(output.read_bytes(), Path(self.env[fixture]).read_bytes())
        self.assertFalse(gene.with_suffix('.pkl').exists() or exon.with_suffix('.pkl').exists())
        self.assertEqual({path: path.read_bytes() for path in protected}, before)
        vcf = self.env['TEST_VCF']
        self.assertEqual(self.priva_calls(), [
            ['clinvar-contract', 'af-paths-check', '--vcf', vcf, '--config', str(config),
             '--gene', str(gene), '--exon', str(exon)],
            ['clinvar-contract', 'header', '--vcf', vcf, '--stage', 'af'],
            ['clinvar-gene-af', 'build', '--vcf', vcf, '--output', str(gene) + '.partial.fixture', '--threads', '1'],
            ['clinvar-contract', 'cache', '--input', str(gene) + '.partial.fixture', '--kind', 'gene-af'],
            ['clinvar-exon-af', 'build', '--vcf', vcf, '--output', str(exon) + '.partial.fixture', '--threads', '1'],
            ['clinvar-contract', 'cache', '--input', str(exon) + '.partial.fixture', '--kind', 'exon-af']])

    def test_af_managed_archives_are_never_rebuilt_in_place(self):
        self.install_vcf('complete')
        managed = {}
        for key in ('TEST_GENE_AF', 'TEST_EXON_AF'):
            managed[key] = self.cache / ('managed.' + key.lower() + '.v1.rkyv')
            managed[key].write_bytes(b'published archive')
            Path(str(managed[key]) + '.source.json').write_text('{"origin": "direct"}\n')
            self.env[key] = str(managed[key])
        result, events = self.run_function('ClinVar_patho_AF_stat config.yaml')
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertEqual(events[2:], [
            'update clinvar_patho_af_stat ' + str(self.cache / 'clinvar.hg19.vep.patho_af_stat.rkyv'),
            'update clinvar_patho_exon_af_stat ' + str(self.cache / 'clinvar.hg19.vep.patho_exon_af_stat.rkyv')])
        self.assertTrue(all(path.read_bytes() == b'published archive' for path in managed.values()))

    def test_af_builder_without_csq_output_fails_before_exon_builder(self):
        self.install_vcf('complete')
        result, events = self.run_function('ClinVar_patho_AF_stat config.yaml', fail='gene_af_no_csq')
        self.assertNotEqual(result.returncode, 0)
        self.assertIn('priva clinvar-gene-af wrote no archive', result.stdout)
        self.assertEqual(events, ['gene_af_builder'])
        self.assertEqual([call[:2] for call in self.priva_calls()], [
            ['clinvar-contract', 'af-paths-check'], ['clinvar-contract', 'header'], ['clinvar-gene-af', 'build']])

    def test_af_stage_rejects_missing_population_annotations(self):
        self.install_vcf('vep')
        result, events = self.run_function('ClinVar_patho_AF_stat config.yaml')
        self.assertNotEqual(result.returncode, 0)
        self.assertEqual(events, [])

    def test_exon_failure_does_not_publish_partial_af_pair(self):
        self.install_vcf('complete')
        result, events = self.run_function('ClinVar_patho_AF_stat config.yaml', fail='exon_af')
        self.assertNotEqual(result.returncode, 0)
        self.assertEqual(events, ['gene_af_builder', 'exon_af_builder'])
        self.assertFalse((self.cache / 'clinvar.hg19.vep.patho_af_stat.rkyv').exists())

    def test_valid_af_outputs_publish_after_original_gene_then_exon_order(self):
        self.install_vcf('complete')
        result, events = self.run_function('ClinVar_patho_AF_stat config.yaml')
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertEqual([event.split()[0] for event in events],
                         ['gene_af_builder', 'exon_af_builder', 'update', 'update'])

    def test_explicit_regeneration_rejects_baseline_output(self):
        before = self.baseline.read_bytes()
        result, events = self.run_function('ClinVar_Gene_stat_regenerate config.yaml summary.tsv mane.tsv "$DATA_DIR/ClinVar/clinvar_2star_stats.pkl.gz"')
        self.assertNotEqual(result.returncode, 0)
        self.assertIn('Candidate output must end in .rkyv', result.stdout)
        self.assertEqual(events, [])
        self.assertEqual(self.priva_calls(), [])
        self.assertEqual(before, self.baseline.read_bytes())

    def regenerate(self, seed=None):
        self.install_vcf('complete')
        for filename in ('summary.tsv', 'mane.tsv'):
            (self.directory / filename).write_text('explicit pinned fixture\n')
        self.env['TEST_DIRECTORY'] = str(self.directory)
        if seed is not None:
            self.env['TEST_SEED'] = seed
        return self.run_function('ClinVar_Gene_stat_regenerate config.yaml "$TEST_DIRECTORY/summary.tsv" "$TEST_DIRECTORY/mane.tsv" "$TEST_DIRECTORY/candidate.rkyv"')

    def test_explicit_regeneration_passes_inputs_and_preserves_baseline(self):
        before = self.baseline.read_bytes()
        result, events = self.regenerate()
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertEqual(len(events), 1)
        self.assertIn('--clinvar-summary ' + str(self.directory / 'summary.tsv'), events[0])
        self.assertIn('--mane-summary ' + str(self.directory / 'mane.tsv'), events[0])
        self.assertEqual(before, self.baseline.read_bytes())
        self.assertTrue((self.directory / 'candidate.rkyv').is_file())
        vcf, temporary = self.env['TEST_VCF'], str(self.directory / 'candidate.partial.fixture.rkyv')
        self.assertEqual(self.priva_calls(), [
            ['clinvar-contract', 'header', '--vcf', vcf, '--stage', 'gene'],
            ['clinvar-gene-stats', 'build', '--vcf', vcf, '--output', temporary, '--python-hash-seed', '0',
             '--clinvar-summary', str(self.directory / 'summary.tsv'),
             '--mane-summary', str(self.directory / 'mane.tsv'), '--threads', '1', '--config', 'config.yaml'],
            ['clinvar-contract', 'cache', '--input', temporary, '--kind', 'gene']])

    def test_explicit_regeneration_uses_the_configured_hash_seed(self):
        result, _ = self.regenerate('17')
        self.assertEqual(result.returncode, 0, result.stderr)
        build = [call for call in self.priva_calls() if call[0] == 'clinvar-gene-stats']
        self.assertEqual(build[0][build[0].index('--python-hash-seed') + 1], '17')
        self.events.write_text('')
        result, events = self.regenerate('-1')
        self.assertNotEqual(result.returncode, 0)
        self.assertIn('priva_python_hash_seed must be a non-negative integer', result.stdout)
        self.assertEqual(events, [])
        self.assertEqual(self.priva_calls(), [])


if __name__ == '__main__':
    unittest.main()
