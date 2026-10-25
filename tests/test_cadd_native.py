"""Portable CADD deployment/routing fixtures; tiny fake assets, no model inference."""
import argparse
import gzip
import importlib.util
import io
import json
import os
from pathlib import Path
import shutil
import subprocess
import tempfile
import unittest
from unittest import mock

ROOT = Path(__file__).resolve().parents[1]
spec = importlib.util.spec_from_file_location('native_cadd_fixture', ROOT/'scripts/cadd_native.py')
N = importlib.util.module_from_spec(spec); spec.loader.exec_module(N)
T, S = N.T, N.S


class NativeCadd(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory(prefix='native_cadd_fixture_', dir=os.environ.get('PRIVA_TEST_TMPDIR', os.environ.get('TMPDIR')))
        self.addCleanup(self.tmp.cleanup); self.root = Path(self.tmp.name)
        self.cache = self.root/'cache'; self.source = self.cache/'CADD/CADD-scripts-1.7.2'
        self.runtime = self.cache/'runtime/cadd-1.7.2'; self.work = self.root/'configured scratch'
        self.work.mkdir(); self.config = dict(assembly='hg38', vep_installed_version='116',
            cache_root=str(self.cache), tmp_dir=str(self.work), cadd_base_dir=str(self.source))
        self.tools = self.root/'tools'; self.tools.mkdir()
        for name in ('bgzip', 'tabix'):
            self.write(self.tools/name, b'#!/bin/sh\n# Tiny fixture: never executed.\nexit 99\n').chmod(0o755)
        self.which = mock.patch('shutil.which', side_effect=lambda name: str(self.tools/name))
        self.which.start(); self.addCleanup(self.which.stop)
        self.seed()

    def write(self, path, data=b'tiny fake asset; not a scientific model\n'):
        path.parent.mkdir(parents=True, exist_ok=True); path.write_bytes(data); return path

    def seed(self):
        annotation = self.source/'data/annotations/GRCh38_v1.7'
        for name in ('regseq/Hyperopt400InclNegatives.json', 'regseq/Hyperopt400InclNegatives.h5',
                     'regseq/reference.fa', 'regseq/reference.fa.genome', 'reference/reference.fa'):
            self.write(annotation/name)
        self.write(self.source/'config/references_GRCh38_v1.7.cfg')
        self.write(self.source/'config/config_GRCh38_v1.7.yml', b'PrescoredFolder: data/prescored/GRCh38_v1.7/incl_anno\n')
        for name in ('src/scripts/VCF2vepVCF.py', 'src/scripts/extract_scored.py', 'Snakefile'): self.write(self.source/name)
        self.prescored = self.source/'data/prescored/GRCh38_v1.7/incl_anno'
        for name in ('b.tsv.gz', 'a.tsv.gz'):
            self.write(self.prescored/name); self.write(self.prescored/(name+'.tbi'))
        hub = self.runtime/'esm-hub-GRCh38'; checkpoints = []
        for i in range(1, 6):
            name = 'esm1v_t33_650M_UR90S_'+str(i)+'.pt'; target = self.write(annotation/'esm'/name)
            link = hub/'checkpoints'/name; link.parent.mkdir(parents=True, exist_ok=True); link.symlink_to(target)
            checkpoints.append(dict(link=str(link), target=str(target)))
        wrapper = self.write(self.runtime/'bin/cadd-native-exec', b'#!/bin/sh\n# Tiny fixture interpreter; never executed.\nexit 99\n')
        wrapper.chmod(0o755); routes = {}
        environments = {}
        for name in ('cadd-python2', 'cadd-esm-python', 'cadd-regseq-python', 'cadd-mmsplice-python'):
            path = wrapper.parent/name; path.symlink_to(wrapper.name); routes[name] = str(path)
            prefix = self.runtime/'envs'/name
            self.write(prefix/'bin/python'); self.write(prefix/'conda-meta/history')
            package = self.write(prefix/'conda-meta/tiny-fixture.json', b'{"fixture_only":true}\n')
            ready = self.runtime/'native-state'/(name+'.ready.json')
            self.write(ready, N.encoded(dict(conda_metadata_sha256={package.name:T.digest(package)})))
            environments[name] = dict(prefix=str(prefix), receipt=str(ready), sha256=T.digest(ready))
        program = self.runtime/'bin/cadd-python-calculate'
        shutil.copy2(ROOT/'scripts/cadd_python_calculate.sh', program); program.chmod(0o755)
        assets = []
        for name in ('config/impute_GRCh38_v1.7.cfg', 'data/models/GRCh38/CADDv1.7-GRCh38.mod',
                     'data/models/GRCh38/conversionTable_CADDv1.7-GRCh38.txt'):
            path = self.write(self.source/name); assets.append(dict(path=str(path), sha256=T.digest(path)))
        peptide = self.write(self.cache/'references/Ensembl116/GRCh38/views/peptide.fa')
        gtf = self.write(peptide.parent/'transcripts.gtf')
        env_candidate = dict(CADD_ESM_TRANSCRIPTS=str(peptide), CADD_MMSPLICE_GTF=str(gtf))
        candidate_path = peptide.parent.parent/'candidate.fixture.ready.json'
        self.write(candidate_path, N.encoded(dict(stage='cadd-selected-transcript-candidate', scientific_qualified=False, cadd_enabled=False)))
        candidate = dict(path=str(candidate_path), sha256=T.digest(candidate_path), assembly='GRCh38',
                         vep_release='116', environment_candidates=env_candidate)
        self.record = dict(backend='native', selection='both', assembly='GRCh38', source_commit=T.PIN,
            status='installed_not_scientifically_qualified', source=str(self.source), fingerprint='a'*64,
            selected_transcript_candidate=candidate,
            runtime=dict(wrappers=dict(path=str(wrapper), sha256=T.digest(wrapper), entry_points=routes), environments=environments),
            assets=dict(esm_hub=dict(path=str(hub), status='ready', checkpoints=checkpoints)),
            calculation_runtime=dict(program=str(program), sha256=T.PYTHON_CALCULATION_SHA256,
                implementation='all_four_calculation_stages_upstream_python', environment=dict(CADD_SOURCE_ROOT=str(self.source),
                CADD_PYTHON2=routes['cadd-python2'], CADD_IMPUTE_CONFIG=assets[0]['path'], CADD_MODEL=assets[1]['path'],
                CADD_PHRED_TABLE=assets[2]['path'])), calculation_assets=assets)
        self.receipt = self.runtime/'receipts/fixture.ready.json'; self.write(self.receipt, N.encoded(self.record))

    def install(self):
        result = N.install_profile(self.config, self.receipt)
        return Path(result['path']), T.read_json(result['path'])

    def test_cold_profile_is_complete_bound_and_unqualified_without_processes(self):
        with mock.patch('subprocess.run', side_effect=AssertionError('No probes/models/downloads')):
            path, profile = self.install()
        self.assertTrue(path.is_relative_to(self.cache))
        self.assertEqual(profile['calculation']['program'], self.record['calculation_runtime']['program'])
        self.assertEqual(profile['calculation']['program_sha256'], '18d10c4533fc9275a12d0ac675b55903b2ee835d4c6f31de2a3a86452c26382d')
        self.assertFalse(profile['production_activation']); self.assertFalse(profile['scientific_qualification'])
        shared = profile['annotation']['shared_vep']
        for key in ('qualification', 'selection_context', 'transcript_policy'): self.assertIsNone(shared[key])
        self.assertIsNone(profile['prescoring']['qualification'])
        self.assertEqual([Path(row['path']).name for row in profile['prescoring']['files']], ['a.tsv.gz', 'b.tsv.gz'])
        self.assertNotIn('CADD_VEP', profile['annotation']['env'])
        self.assertNotIn('standalone', profile['calculation'])
        self.assertEqual(N.check(self.config, path), profile)

    def test_warm_install_and_verify_reuse_exact_bytes_without_network_or_processes(self):
        path, profile = self.install(); directory = path.parent
        before = {p.name:(p.stat().st_ino, p.stat().st_mtime_ns, p.read_bytes()) for p in directory.iterdir()}
        with mock.patch('subprocess.run', side_effect=AssertionError('No subprocess')), mock.patch('urllib.request.urlopen', side_effect=AssertionError('No network')):
            self.assertEqual(N.install_profile(self.config, self.receipt)['path'], str(path))
            N.install_profile(self.config, self.receipt, verify_only=True)
        self.assertEqual(before, {p.name:(p.stat().st_ino, p.stat().st_mtime_ns, p.read_bytes()) for p in directory.iterdir()})

    def test_unqualified_execution_fails_before_any_output_or_scientific_process(self):
        path, profile = self.install(); output = self.root/'must stay absent'
        with mock.patch.object(T, 'execute') as execute:
            with self.assertRaisesRegex(S.SharedVepError, 'unqualified'):
                N.run(self.config, path, self.root/'unprovided bundle', output)
        execute.assert_not_called(); self.assertFalse(output.exists()); self.assertEqual(list(self.work.iterdir()), [])

    def test_missing_deployment_verify_does_not_publish(self):
        with self.assertRaisesRegex(T.ContractError, 'incomplete'): N.install_profile(self.config, self.receipt, verify_only=True)
        self.assertFalse((self.runtime/'profiles').exists())

    def test_corrupt_existing_profile_retained_and_refused(self):
        path, profile = self.install(); path.write_text('unrelated existing contents\n')
        with self.assertRaisesRegex(T.ContractError, 'differs'): N.install_profile(self.config, self.receipt)
        self.assertEqual(path.read_text(), 'unrelated existing contents\n')

    def test_changed_alias_checkpoint_program_and_prescore_membership_rejected(self):
        path, profile = self.install()
        checkpoint = Path(self.record['assets']['esm_hub']['checkpoints'][0]['link']); target = checkpoint.resolve()
        checkpoint.unlink(); checkpoint.symlink_to(self.record['assets']['esm_hub']['checkpoints'][1]['target'])
        with self.assertRaisesRegex(T.ContractError, 'checkpoint route'): N.check(self.config, path)
        checkpoint.unlink(); checkpoint.symlink_to(target)
        program = Path(self.record['calculation_runtime']['program']); original = program.read_bytes(); program.write_bytes(b'changed')
        with self.assertRaisesRegex(T.ContractError, 'calculation binding'): N.check(self.config, path)
        program.write_bytes(original)
        self.write(self.prescored/'late.tsv.gz'); self.write(self.prescored/'late.tsv.gz.tbi')
        with self.assertRaisesRegex(T.ContractError, 'membership|prescore order'): N.check(self.config, path)

    def test_native_ready_receipt_and_conda_metadata_drift_are_rejected(self):
        path, profile = self.install()
        environment = self.record['runtime']['environments']['cadd-esm-python']
        ready = Path(environment['receipt']); original = ready.read_bytes(); ready.write_bytes(b'changed receipt')
        with self.assertRaisesRegex(T.ContractError, 'environment receipt'): N.check(self.config, path)
        ready.write_bytes(original)
        metadata = Path(environment['prefix'])/'conda-meta/tiny-fixture.json'
        original = metadata.read_bytes(); metadata.write_bytes(b'changed package')
        with self.assertRaisesRegex(RuntimeError, 'package metadata'): N.check(self.config, path)
        metadata.write_bytes(original)
        self.write(metadata.parent/'unexpected-package.json', b'{}\n')
        with self.assertRaisesRegex(RuntimeError, 'package membership'): N.check(self.config, path)

    def test_native_source_identity_has_no_nextflow_dependency_or_implicit_fallback(self):
        identity = T.source_identity('native')
        self.assertIn('scripts/cadd_native.py', identity)
        self.assertIn('scripts/cadd_runtime.py', identity)
        self.assertFalse(any(name.endswith('.nf') or name=='nextflow.config' for name in identity))
        self.assertIn('bin/cadd/CaddSelection.pm', identity)
        self.assertIn('bin/cadd/upstream_calculate.sh', identity)
        with self.assertRaisesRegex(T.ContractError, 'scope'): T.source_identity('guess')

    def test_imported_runtime_helper_digest_drift_fails_preflight_before_stages(self):
        path, profile = self.install(); digest = T.digest
        helper = T.ROOT/'scripts/cadd_runtime.py'
        def changed_digest(value):
            return '0'*64 if Path(value) == helper else digest(value)
        # Isolate the engineering source binding; these fake assets are never
        # qualified or executed, and actual source files remain untouched.
        with mock.patch.object(S, 'verify_profile'), mock.patch.object(N, 'verify_prescoring_profile'), \
                mock.patch.object(T, 'digest', side_effect=changed_digest), mock.patch.object(T, 'stage_preflight') as stages:
            with self.assertRaisesRegex(T.ContractError, 'adapter source changed'):
                N.preflight(self.config, path)
        stages.assert_not_called(); self.assertEqual(list(self.work.iterdir()), [])

    def test_profile_cannot_check_unused_interpreters_while_executing_others(self):
        path, profile = self.install()
        for key in ('CADD_PYTHON2', 'CADD_ESM_PYTHON', 'CADD_REGSEQ_PYTHON', 'CADD_MMSPLICE_PYTHON'):
            changed = json.loads(json.dumps(profile)); changed['annotation']['env'][key] = str(self.tools/'bgzip')
            changed_path = self.root/(key+'.json'); changed_path.write_bytes(N.encoded(changed))
            with self.subTest(key=key), self.assertRaisesRegex(T.ContractError, 'interpreter binding'):
                N.check(self.config, changed_path)

    def test_installed_original_program_is_selected_and_corruption_rejected(self):
        path, profile = self.install()
        with mock.patch.object(T.subprocess, 'run') as command:
            T.run_backend(profile, 'calculation', self.root/'annotations.tsv.gz', self.root/'scores.tsv', self.work)
        self.assertEqual(command.call_args.args[0][1], self.record['calculation_runtime']['program'])
        self.assertEqual(command.call_args.kwargs['env']['CADD_PYTHON2'], self.record['calculation_runtime']['environment']['CADD_PYTHON2'])
        with mock.patch.object(T, 'verify_source'), mock.patch.object(T, 'verify_resources', return_value={}):
            T.stage_preflight(profile, 'calculation')
            Path(profile['calculation']['program']).write_bytes(b'suspended fake backend')
            with self.assertRaisesRegex(T.ContractError, 'pinned original Python'): T.stage_preflight(profile, 'calculation')

    def selected_fixture(self):
        cols = S.REQUIRED_CSQ; values = {key:'' for key in cols}
        values.update(Allele='G', Consequence='missense_variant', Feature_type='Transcript', Gene='G1', Feature='T1', SYMBOL='G1')
        vcf = self.root/'fake-selected.vcf'
        self.write(vcf, ('##fileformat=VCFv4.2\n##INFO=<ID=CSQ,Number=.,Type=String,Description="Format: '+ '|'.join(cols)+'">\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n17\t10\t.\tA\tG\t.\t.\tCSQ='+ '|'.join(values[x] for x in cols)+'\n').encode())
        return vcf

    def split_fixture(self, path, profile, selected):
        """Closed synthetic proof for routing tests only, not qualification evidence."""
        evidence = self.write(self.root/'fake-policy-evidence.json', b'fixture only\n')
        qualification = self.root/'fake-policy-qualification.json'
        self.write(qualification, N.encoded(dict(schema_version=1, status='qualified', policy_id=N.PRESCORE_POLICY,
            profile_id=profile['profile_id'], resource_manifest_sha256=profile['contract']['resource_manifest_sha256'],
            stub=False, checks=N.PRESCORE_CHECKS, evidence=N.file_pin(evidence))))
        profile['prescoring']['qualification'] = N.file_pin(qualification); path.write_bytes(N.encoded(profile))
        mapping = self.write(self.root/'source-occurrences.json', b'{"fixture_only":true}\n')
        prepared = self.root/'prepared.vcf'
        with selected.open() as incoming, prepared.open('w') as outgoing:
            outgoing.write('#CHROM\tPOS\tID\tREF\tALT\n')
            for line in incoming:
                if not line.startswith('#'): outgoing.write('\t'.join(line.rstrip('\n').split('\t')[:5])+'\n')
        previous = N.file_pin(prepared); steps = []
        for i, item in enumerate(profile['prescoring']['files']):
            found = self.write(self.root/('found%d.tsv'%i), b'#Fake publisher header\n')
            remaining = self.write(self.root/('remaining%d.vcf'%i), prepared.read_bytes())
            step = dict(publisher_path=item['path'], input=previous, found=N.file_pin(found), remaining=N.file_pin(remaining))
            steps.append(step); previous=step['remaining']
        key_path = self.root/'novel-keys.tsv'; count = 0
        with key_path.open('w') as stream:
            for key in S.variant_keys(selected): stream.write('\t'.join(key)+'\n'); count += 1
        keys = dict(**N.file_pin(key_path), count=count)
        proof = dict(schema_version=1, status='completed', stage='upstream-prescore-novel-subset-v1',
            policy_id=N.PRESCORE_POLICY, profile_sha256=T.digest(path), prescore_inputs_sha256=N.json_hash(profile['prescoring']),
            full_source=N.file_pin(selected), mapping=N.file_pin(mapping), prepared=N.file_pin(prepared),
            steps=steps, novel=previous, novel_keys=keys)
        proof_path = self.write(self.root/'publisher-split.json', N.encoded(proof))
        bundle = self.root/'fake bundle'; bundle.mkdir()
        with gzip.open(bundle/'vep.vcf.gz','wb') as stream: stream.write(selected.read_bytes())
        pipeline = {key:proof[key] for key in ('stage','profile_sha256','prescore_inputs_sha256','full_source','mapping')}
        pipeline.update(prescore_split=N.file_pin(proof_path), novel_keys_sha256=keys['sha256'])
        self.write(bundle/'manifest.json', N.encoded(dict(pipeline_provenance=pipeline)))
        return bundle, proof_path, proof

    def test_generic_selected_bundle_cannot_enter_true_novel_native_endpoint(self):
        path, profile = self.install(); selected = self.selected_fixture()
        bundle, proof_path, proof = self.split_fixture(path, profile, selected)
        (bundle/'manifest.json').write_text('{}\n'); output = self.root/'must stay absent'
        with mock.patch.object(N,'preflight',return_value=profile), mock.patch.object(S,'verify',return_value=(selected,{})), mock.patch.object(T,'execute') as execute:
            with self.assertRaisesRegex(T.ContractError, 'profile-bound original publisher-prescore split proof'):
                N.run(self.config,path,bundle,output)
        execute.assert_not_called(); self.assertFalse(output.exists()); self.assertEqual(list(self.work.iterdir()),[])

    def test_split_proof_requires_all_files_order_and_exact_novel_keys(self):
        path, profile = self.install(); selected = self.selected_fixture()
        bundle, proof_path, proof = self.split_fixture(path, profile, selected)
        self.assertEqual(len(Path(proof['novel']['path']).read_text().splitlines()[-1].split('\t')),5)
        N.verify_novel_input(profile,path,bundle)
        for change, message in (('omit','omitted'),('order','order'),('keys','true-novel')):
            broken = json.loads(json.dumps(proof))
            if change=='omit': broken['steps'].pop()
            elif change=='order': broken['steps'][0]['publisher_path']=profile['prescoring']['files'][1]['path']
            else:
                key_path=Path(broken['novel_keys']['path']); key_path.write_text('17\t999\tA\tG\n')
                broken['novel_keys']['sha256']=T.digest(key_path)
            proof_path.write_bytes(N.encoded(broken))
            receipt=T.read_json(bundle/'manifest.json'); receipt['pipeline_provenance']['prescore_split']=N.file_pin(proof_path)
            receipt['pipeline_provenance']['novel_keys_sha256']=broken['novel_keys']['sha256']
            (bundle/'manifest.json').write_bytes(N.encoded(receipt))
            with self.subTest(change=change), self.assertRaisesRegex(T.ContractError,message): N.verify_novel_input(profile,path,bundle)

    def test_streamed_key_proof_preserves_large_fixture_order_and_occurrences(self):
        path, profile = self.install(); selected = self.selected_fixture()
        text=selected.read_text().splitlines(); template=text[-1].split('\t')
        with selected.open('w') as stream:
            stream.write('\n'.join(text[:-1])+'\n')
            for i in range(20000):
                row=template[:]; row[1]=str(10+i//2); stream.write('\t'.join(row)+'\n')
        bundle, proof_path, proof = self.split_fixture(path,profile,selected)
        self.assertEqual(proof['novel_keys']['count'],20000)
        self.assertIsInstance(proof['novel_keys'],dict)
        N.verify_novel_input(profile,path,bundle)

    def test_native_composition_preserves_complete_rows_missingness_and_configured_scratch(self):
        path, profile = self.install(); selected = self.selected_fixture(); calls = []
        columns = profile['contract']['columns']; row = ['NA']*len(columns)
        for name,value in (('Chrom','17'),('Pos','10'),('Ref','A'),('Alt','G'),('FeatureID','T1')): row[columns.index(name)] = value
        def backend(p, stage, source, output, work, annotation_source='shared-vep'):
            calls.append((stage, Path(work), annotation_source))
            self.assertTrue(Path(work).is_relative_to(self.work))
            if stage=='annotation':
                with gzip.open(output,'wt') as stream: stream.write('#'+'\t'.join(columns)+'\n'+'\t'.join(row)+'\n')
            else:
                self.assertEqual(p['calculation']['program'], self.record['calculation_runtime']['program'])
                Path(output).write_text('#'+'\t'.join(columns+['RawScore','PHRED'])+'\n'+'\t'.join(row+['1.000000','10'])+'\n')
        bundle, _, _ = self.split_fixture(path, profile, selected); output = self.root/'result'
        # Only this routing fixture bypasses qualification/package checks. No
        # fake asset is represented as scientifically qualified or executed.
        with mock.patch.object(N,'preflight',return_value=profile), mock.patch.object(S,'verify',return_value=(selected,{})), \
                mock.patch.object(T,'stage_preflight'), mock.patch.object(T,'run_backend',side_effect=backend):
            result = N.run(self.config,path,bundle,output)
        self.assertEqual([item[0] for item in calls], ['annotation','calculation'])
        self.assertEqual(result['execution_status'],'completed'); self.assertFalse(result['scientific_qualification'])
        self.assertTrue((output/'annotation_bundle/manifest.json').is_file())
        self.assertTrue((output/'calculation_bundle/manifest.json').is_file())
        score = (output/'calculation_bundle/scores.tsv').read_text().splitlines()[1].split('\t')
        self.assertEqual(score[:-2],row); self.assertEqual(len(score),153)
        self.assertFalse((output/'annotation_bundle.work').exists())

    def test_process_failure_exit_is_preserved_and_no_completion_manifest_is_published(self):
        path, profile = self.install(); selected = self.selected_fixture()
        bundle, _, _ = self.split_fixture(path,profile,selected); config_path = self.root/'config.json'
        config_path.write_text(json.dumps(self.config)); output = self.root/'failed result'
        with mock.patch.object(N,'preflight',return_value=profile), mock.patch.object(S,'verify',return_value=(selected,{})), \
                mock.patch.object(T,'stage_preflight'), mock.patch.object(T,'run_backend',side_effect=subprocess.CalledProcessError(37,['tiny fake backend'])), \
                mock.patch('sys.stderr',new_callable=io.StringIO):
            status=N.main(['run','--config',str(config_path),'--profile',str(path),'--shared-vep-bundle',str(bundle),'--output',str(output)])
        self.assertEqual(status,37)
        self.assertEqual(T.read_json(output/'run.json')['execution_status'],'failed')
        self.assertEqual(T.read_json(output/'run.json')['process_exit_code'],37)
        self.assertFalse((output/'annotation_bundle/manifest.json').exists())
        self.assertFalse((output/'calculation_bundle').exists())


if __name__ == '__main__': unittest.main()
