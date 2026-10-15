#!/usr/bin/env python3
"""Mocked offline integration: no biological inference or cache downloads."""
import copy
import contextlib
import configparser
import io
import gzip
import hashlib
import importlib.util
import json
import os
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest
from unittest.mock import patch

ROOT = Path(__file__).resolve().parents[2]
SPEC = importlib.util.spec_from_file_location('cadd_task', ROOT / 'bin/cadd_task.py')
C = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(C)

MOCK = r'''
from pathlib import Path
import gzip,json,os,sys
args=sys.argv[1:]
name=Path(args[0]).name if args and args[0].endswith('.py') else Path(sys.argv[0]).name
if name=='git':
    if args[2:]==['rev-parse','HEAD']:
        print(os.environ.get('MOCK_REVISION','1fead0c1e02d3f0a8599bdb04e6b3a2dd70b4866'));sys.exit(0)
    if args[2:]==['diff','--quiet','HEAD','--','src','config']:
        sys.exit(int(os.environ.get('MOCK_DIRTY','0')))
    if args[2:]==['ls-files','--others','--exclude-standard','--','src','config']:
        print(os.environ.get('MOCK_UNTRACKED',''));sys.exit(0)
    sys.exit(129)
def flag(key):return args[args.index(key)+1]
def log():
    with open(os.environ['MOCK_EVENTS'],'a') as f:f.write(name+'\n')
if name==os.environ.get('MOCK_FAIL_STAGE'):sys.exit(23)
if name in ['max_line_hierarchy.py','appendPHREDscore.py']:
    contents=sys.stdin.read()
    log()
    if name=='max_line_hierarchy.py':
        assert '--all' in args
        sys.stdout.write(contents)
    else:
        assert Path(flag('-t')).is_file()
        for line in contents.splitlines():
            print(line+'\t'+('PHRED' if line.startswith('#') else '12.3'))
    sys.exit(0)
log()
if name in ['VCF2vepVCF.py','vep']:sys.stdout.write(sys.stdin.read())
elif name=='annotateVEPvcf.py':
    sys.stdin.read();print('##metadata')
    columns=json.loads(Path(os.environ['MOCK_COLUMNS']).read_text())
    print('#'+'\t'.join(columns))
    if not os.environ.get('MOCK_EMPTY'):
        fields=['NA']*len(columns)
        for key,value in [('Chrom','1'),('Pos','100'),('Ref','A'),('Alt','G'),('Type','SNV'),('Consequence','NON_SYNONYMOUS'),('FeatureID',os.environ.get('MOCK_FEATURE_ID','NA'))]:fields[columns.index(key)]=value
        print('\t'.join(fields))
elif name=='trackTransformation.py':
    text=sys.stdin.read()
    assert text.startswith('#Chrom\t') and '##metadata' not in text
    with gzip.open(flag('-o'),'wt') as f:
        for line in text.splitlines()[1:]:f.write('1,2,3\n')
elif name=='predictSKmodel.py':
    with gzip.open(flag('-a'),'rt') as f:
        for index,line in enumerate(f):
            if os.environ.get('MOCK_DROP_ROW') and index:continue
            print(line.rstrip()+'\t'+('RawScore' if index==0 else '1.500000'))
elif name=='MMSplice.py':
    with gzip.open(flag('-i'),'rt') as f:sys.stdout.write(f.read())
elif name=='tabix':Path(args[-1]+'.tbi').touch()
else:
    source=flag('--variants') if '--variants' in args else flag('--input')
    Path(flag('--output')).write_bytes(Path(source).read_bytes())
'''


class UpstreamCadd(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory(prefix='cadd-tests-', dir=os.environ.get('TMPDIR'))
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)
        self.profile = self.root / 'profile.json'
        self.p = json.loads((ROOT / 'conf/cadd/GRCh37.reference.json').read_text())
        self.p['software_id'] = 'mock-only-not-scientific'
        self.bin = self.root / 'bin'
        self.bin.mkdir()
        for name in ['python', 'git', 'vep', 'tabix']:
            path = self.bin / name
            path.write_text('#!' + sys.executable + '\n' + MOCK)
            path.chmod(0o755)
        self.source = self.root / 'source'
        for name in ['VCF2vepVCF.py', 'annotateVEPvcf.py', 'trackTransformation.py', 'predictSKmodel.py',
                     'max_line_hierarchy.py', 'appendPHREDscore.py', 'lib/tools/MMSplice.py',
                     'lib/tools/esmScore/esmScore_missense_av_fast.py', 'lib/tools/esmScore/esmScore_frameshift_av.py',
                     'lib/tools/esmScore/esmScore_inFrame_av.py', 'lib/tools/regulatorySequence/predictVariants.py']:
            path = self.source / 'src/scripts' / name
            path.parent.mkdir(parents=True, exist_ok=True)
            path.write_text('# Mock placeholder, never executed.\n')
        self.assets = self.root / 'assets'
        self.assets.mkdir()
        for name in ['peptides.fa', 'regseq.json', 'weights.h5', 'reference.fa', 'reference.fa.fai',
                     'genome', 'references.cfg', 'transcripts.gtf', 'impute.cfg', 'model.mod', 'phred.tsv']:
            (self.assets / name).write_text('Mock asset\n')
        (self.assets / 'vep').mkdir()
        (self.assets / 'vep/cache').write_text('Mock cache\n')
        (self.assets / 'esm/checkpoints').mkdir(parents=True)
        for number in range(1, 6):
            (self.assets / f'esm/checkpoints/esm1v_t33_650M_UR90S_{number}.pt').write_text('Mock weights\n')
        env = self.p['annotation']['env']
        for key in ['CADD_PYTHON2', 'CADD_ESM_PYTHON', 'CADD_REGSEQ_PYTHON', 'CADD_MMSPLICE_PYTHON']:
            env[key] = str(self.bin / 'python')
        env.update(CADD_SOURCE_ROOT=str(self.source), CADD_VEP=str(self.bin / 'vep'),
                   CADD_VEP_CACHE=str(self.assets / 'vep'), CADD_ESM_MODEL_DIR=str(self.assets / 'esm'),
                   CADD_ESM_TRANSCRIPTS=str(self.assets / 'peptides.fa'), CADD_REGSEQ_MODEL=str(self.assets / 'regseq.json'),
                   CADD_REGSEQ_WEIGHTS=str(self.assets / 'weights.h5'), CADD_REGSEQ_REFERENCE=str(self.assets / 'reference.fa'),
                   CADD_REGSEQ_GENOME=str(self.assets / 'genome'), CADD_REFERENCE_CONFIG=str(self.assets / 'references.cfg'),
                   CADD_MMSPLICE_GTF=str(self.assets / 'transcripts.gtf'), CADD_REFERENCE_FASTA=str(self.assets / 'reference.fa'),
                   CADD_BGZIP='/usr/bin/gzip', CADD_TABIX=str(self.bin / 'tabix'))
        calc = self.p['calculation']
        calc['env'] = {'CADD_SOURCE_ROOT':str(self.source), 'CADD_PYTHON2':str(self.bin / 'python')}
        for key, name in [('CADD_IMPUTE_CONFIG','impute.cfg'), ('CADD_MODEL','model.mod'), ('CADD_PHRED_TABLE','phred.tsv')]:
            calc['env'][key] = str(self.assets / name)
        calc['verify_files'] = [{'path':str(self.assets / name), 'sha256':C.digest(self.assets / name)}
                                for name in ['impute.cfg','model.mod','phred.tsv']]
        calc['resource_roots'] = [item['path'] for item in calc['verify_files']]
        self.p['annotation']['resource_roots'] = [str(self.assets)]
        self.p['resource_manifest_path'] = str(self.root / 'resources.json')
        self.columns = self.root / 'columns.json'
        self.columns.write_text(json.dumps(self.p['contract']['columns']))
        self.events = self.root / 'events'
        self.environment = patch.dict(os.environ, {'PATH':str(self.bin)+os.pathsep+os.environ['PATH'],
            'MOCK_EVENTS':str(self.events), 'MOCK_COLUMNS':str(self.columns)}, clear=False)
        self.environment.start()
        self.addCleanup(self.environment.stop)
        self.vcf = self.root / 'input.vcf'
        self.vcf.write_text('##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n1\t100\t.\tA\tG\t.\t.\t.\n')
        self.make_reference_config()
        self.save()
        self.lock()

    def make_reference_config(self):
        self.p['annotation']['required_annotations'] = json.loads((ROOT/'conf/cadd'/(self.p['contract']['assembly']+'.reference.json')).read_text())['annotation']['required_annotations']
        cfg = configparser.ConfigParser()
        cfg['Path'] = {'roottype':'path','rootdir':str(self.assets),'annodir':'','reference':'/reference.fa'}
        cfg['Annotations'] = {}
        for index,(name,kind) in enumerate(self.p['annotation']['required_annotations'].items()):
            if kind == 'computed':
                cfg['Annotations'][name] = 'True'
            else:
                path = self.assets / ('track-' + str(index) + ('.gz' if kind=='tabix' else '.txt'))
                path.write_text('Mock track\n' if kind!='tabix-list' else 'mock.list.member.gz\n')
                cfg['Annotations'][name] = '/' + path.name
                if kind=='tabix':Path(str(path)+'.tbi').write_text('Mock index\n')
                if kind=='tabix-list':
                    (self.assets/'mock.list.member.gz').write_text('Mock track\n')
                    (self.assets/'mock.list.member.gz.tbi').write_text('Mock index\n')
        with open(self.assets/'references.cfg','w') as handle:cfg.write(handle)

    def save(self):
        self.profile.write_text(json.dumps(self.p))

    def lock(self):
        target = Path(self.p['resource_manifest_path'])
        if target.exists():target.unlink()
        with contextlib.redirect_stdout(io.StringIO()):
            C.lock_resources(argparse_namespace(profile=str(self.profile),output=str(target),checksum=False))
        self.p['contract']['resource_manifest_sha256']=C.digest(target)
        self.save()

    def run_task(self, command, input_path, name, stub=False):
        output=self.root/name
        C.execute(argparse_namespace(command=command,profile=str(self.profile),assembly=self.p['contract']['assembly'],
            input=str(input_path),output=str(output),stub=stub,annotation_source='reference-vcf'))
        return output

    def calls(self):
        return self.events.read_text().splitlines() if self.events.exists() else []

    def test_both_assembly_annotation_branches_and_all_python_calculation(self):
        for assembly in ['GRCh37','GRCh38']:
            with self.subTest(assembly=assembly):
                self.p['contract'].update(assembly=assembly,columns=json.loads((ROOT/'conf/cadd'/f'{assembly}.reference.json').read_text())['contract']['columns'])
                self.columns.write_text(json.dumps(self.p['contract']['columns']))
                self.make_reference_config();self.save();self.lock();self.events.write_text('')
                bundle=self.run_task('annotate',self.vcf,assembly+'-anno')
                expected=['VCF2vepVCF.py','vep','esmScore_missense_av_fast.py','esmScore_frameshift_av.py','esmScore_inFrame_av.py','predictVariants.py']
                if assembly=='GRCh38':expected+=['tabix','MMSplice.py']
                expected+=['annotateVEPvcf.py']
                self.assertEqual(self.calls(),expected)
                with gzip.open(bundle/'annotations.tsv.gz','rt') as f:self.assertFalse(f.read().startswith('##'))
                self.events.write_text('')
                result=self.run_task('calculate',bundle,assembly+'-calc')
                self.assertEqual(self.calls(),['trackTransformation.py','predictSKmodel.py','max_line_hierarchy.py','appendPHREDscore.py'])
                self.assertEqual(json.loads((result/'manifest.json').read_text())['row_count'],1)
                self.assertTrue((result/'scores.tsv').read_text().endswith('\t1.500000\t12.3\n'))
                self.assertFalse((result/'work').exists())

    def test_calculation_without_annotation_resources_or_environment(self):
        bundle=self.run_task('annotate',self.vcf,'anno')
        del self.p['annotation']
        self.save()
        # Keep the shared manifest identity: calculate must never inspect absent annotation paths.
        (self.assets/'peptides.fa').unlink()
        self.run_task('calculate',bundle,'calc')

    def test_source_revision_dirty_or_untracked_are_rejected(self):
        for key,value in [('MOCK_REVISION','wrong'),('MOCK_DIRTY','1'),('MOCK_UNTRACKED','src/rogue.py')]:
            with self.subTest(key=key),patch.dict(os.environ,{key:value}):
                with self.assertRaises((C.ContractError,subprocess.CalledProcessError)):
                    self.run_task('annotate',self.vcf,key)

    def test_model_and_resource_changes_rejected(self):
        (self.assets/'model.mod').write_text('Changed model\n')
        with self.assertRaises(C.ContractError):C.stage_preflight(self.p,'calculation')
        self.lock()
        with self.assertRaisesRegex(C.ContractError,'checksum'):C.stage_preflight(self.p,'calculation')

    def test_incomplete_columns_and_assembly_mismatch_rejected(self):
        with self.assertRaisesRegex(C.ContractError,'assembly'):C.load_profile(self.profile,'hg38')
        self.p['contract']['columns']=self.p['contract']['columns'][:-1];self.save()
        with self.assertRaisesRegex(C.ContractError,'Complete ordered'):C.load_profile(self.profile)

    def test_bundle_contract_and_content_change_rejected(self):
        bundle=self.run_task('annotate',self.vcf,'anno')
        changed=copy.deepcopy(self.p);changed['contract']['resource_manifest_sha256']='a'*64
        with self.assertRaisesRegex(C.ContractError,'mismatch'):C.verify_bundle(bundle,changed)
        with open(bundle/'annotations.tsv.gz','ab') as f:f.write(b'changed')
        with self.assertRaisesRegex(C.ContractError,'content changed'):C.verify_bundle(bundle,self.p)

    def test_failed_science_does_not_publish_completion(self):
        bundle=self.run_task('annotate',self.vcf,'anno')
        with patch.dict(os.environ,{'MOCK_FAIL_STAGE':'predictSKmodel.py'}):
            with self.assertRaises(subprocess.CalledProcessError):self.run_task('calculate',bundle,'failed')
        self.assertFalse((self.root/'failed/manifest.json').exists())

    def test_annotation_row_loss_rejected(self):
        bundle=self.run_task('annotate',self.vcf,'anno')
        with patch.dict(os.environ,{'MOCK_DROP_ROW':'1'}):
            with self.assertRaisesRegex(C.ContractError,'dropped'):self.run_task('calculate',bundle,'dropped')
        self.assertFalse((self.root/'dropped/manifest.json').exists())

    def test_empty_tables_and_import_preserve_contract(self):
        with patch.dict(os.environ,{'MOCK_EMPTY':'1'}):bundle=self.run_task('annotate',self.vcf,'anno')
        imported=self.run_task('import-annotations',bundle/'annotations.tsv.gz','imported')
        result=self.run_task('calculate',imported,'calc')
        self.assertEqual(json.loads((result/'manifest.json').read_text())['row_count'],0)

    def test_stub_artifact_cannot_enter_real_calculation(self):
        self.p['stub']=True;self.save()
        bundle=self.run_task('annotate',self.vcf,'stub',True)
        self.p.pop('stub');self.save()
        with self.assertRaisesRegex(C.ContractError,'Stub'):self.run_task('calculate',bundle,'calc')

    def shared_fixture(self):
        shared = self.root/'shared-vep'
        shared.mkdir()
        cols = C.cadd_shared_vep.REQUIRED_CSQ
        values = {key:'' for key in cols}
        values.update(Allele='G',Consequence='missense_variant',Feature_type='Transcript',Feature='ENST0001')
        text = '##fileformat=VCFv4.2\n##INFO=<ID=CSQ,Number=.,Type=String,Description="Consequence. Format: ' + '|'.join(cols) + '">\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n1\t100\t.\tA\tG\t.\t.\tCSQ=' + '|'.join(values[key] for key in cols) + '\n'
        with gzip.open(shared/'vep.vcf.gz','wt') as handle:handle.write(text)
        original = self.root/'single-priva-vep.vcf'
        original.write_text(text)
        proof = self.root/'selection-proof.json'
        proof.write_text(json.dumps({'scope':'synthetic-test-only','cases':1}))
        qualification = self.root/'selection-qualification.json'
        qualification.write_text(json.dumps({'schema_version':1,'status':'qualified','assembly':'GRCh37',
            'vep_release':'110','transcript_release':'110','selection':'cadd-per-gene','adapter_id':'mock-adapter',
            'stub':False,'checks':{'per_gene_selection':'exact','feature_id_joins':'exact','variant_representation':'exact','missingness':'preserved'},
            'evidence':{'path':str(proof),'sha256':C.digest(proof)}}))
        settings = self.p['annotation']['shared_vep']
        settings.update(adapter_id='mock-adapter',qualification={'path':str(qualification),'sha256':C.digest(qualification)})
        receipt = {'schema_version':1,'stage':'cadd-shared-vep-view','assembly':'GRCh37',
            'profile_id':self.p['profile_id'],'resource_manifest_sha256':self.p['contract']['resource_manifest_sha256'],
            'vep_release':'110','transcript_release':'110','selection':'cadd-per-gene','adapter_id':'mock-adapter',
            'qualification_sha256':C.digest(qualification),'stub':False,
            'source_vep':{'path':str(original),'sha256':C.digest(original)},
            'output_sha256':C.digest(shared/'vep.vcf.gz'),'expected_annotation_rows':1}
        (shared/'manifest.json').write_text(json.dumps(receipt))
        self.save()
        return shared

    def test_shared_vep_never_runs_prepare_or_vep(self):
        shared=self.shared_fixture()
        with patch.dict(os.environ,{'MOCK_FEATURE_ID':'ENST0001'}):
            C.execute(argparse_namespace(command='annotate',profile=str(self.profile),assembly='GRCh37',
                input=str(shared),output=str(self.root/'shared-result'),stub=False,annotation_source='shared-vep'))
        self.assertEqual(self.calls(),['esmScore_missense_av_fast.py','esmScore_frameshift_av.py','esmScore_inFrame_av.py','predictVariants.py','annotateVEPvcf.py'])
        receipt=json.loads((self.root/'shared-result/manifest.json').read_text())
        self.assertEqual(receipt['input_manifest_sha256'],C.digest(shared/'manifest.json'))

    def test_unqualified_shared_vep_and_release_mismatch_rejected(self):
        with self.assertRaisesRegex(C.cadd_shared_vep.SharedVepError,'unqualified'):
            C.cadd_shared_vep.verify(self.root,self.p)
        shared=self.shared_fixture()
        self.p['annotation']['shared_vep']['esm_transcript_release']='113'
        with self.assertRaisesRegex(C.cadd_shared_vep.SharedVepError,'releases differ'):
            C.cadd_shared_vep.verify(shared,self.p)

    def test_shared_vep_feature_join_mismatch_fails_completion(self):
        shared=self.shared_fixture()
        with patch.dict(os.environ,{'MOCK_FEATURE_ID':'WRONG'}):
            with self.assertRaisesRegex(C.cadd_shared_vep.SharedVepError,'FeatureID join'):
                C.execute(argparse_namespace(command='annotate',profile=str(self.profile),assembly='GRCh37',
                    input=str(shared),output=str(self.root/'bad-feature'),stub=False,annotation_source='shared-vep'))
        self.assertFalse((self.root/'bad-feature/manifest.json').exists())

    def test_missing_csq_field_and_changed_source_are_rejected(self):
        shared=self.shared_fixture()
        text=gzip.open(shared/'vep.vcf.gz','rt').read().replace('SIFT','REMOVED')
        with gzip.open(shared/'vep.vcf.gz','wt') as handle:handle.write(text)
        with self.assertRaisesRegex(C.cadd_shared_vep.SharedVepError,'required CSQ'):
            list(C.cadd_shared_vep.expected_features(shared/'vep.vcf.gz'))
        with self.assertRaisesRegex(C.cadd_shared_vep.SharedVepError,'view changed'):
            C.cadd_shared_vep.verify(shared,self.p)

    def test_disabled_source_or_missing_index_is_not_variant_missingness(self):
        cfg=configparser.ConfigParser();cfg.read(self.assets/'references.cfg')
        cfg.set('Annotations','grantham','Ignore')
        with open(self.assets/'references.cfg','w') as handle:cfg.write(handle)
        self.lock()
        with self.assertRaisesRegex(C.ContractError,'resource is disabled'):
            C.stage_preflight(self.p,'annotation')
        self.make_reference_config();self.lock()
        missing=next(self.assets.glob('track-*.gz.tbi'));missing.unlink();self.lock()
        with self.assertRaises(FileNotFoundError):C.stage_preflight(self.p,'annotation')

    def test_import_calculate_profile_needs_no_annotation_resources(self):
        bundle=self.run_task('annotate',self.vcf,'anno')
        del self.p['annotation'];self.save()
        (self.assets/'references.cfg').unlink()
        imported=self.run_task('import-annotations',bundle/'annotations.tsv.gz','imported')
        self.run_task('calculate',imported,'calc')

    def test_output_identity_guard_is_uncached(self):
        args=argparse_namespace(profile=str(self.profile),assembly='GRCh37',stub=False,mode='annotate',
            input=str(self.vcf),output=str(self.root/'published'),annotation_source='reference-vcf')
        with contextlib.redirect_stdout(io.StringIO()):
            C.preflight(args);C.preflight(args)
        self.p['software_id']='different';self.save()
        with self.assertRaisesRegex(C.ContractError,'different inputs/profile/source'):C.preflight(args)


def argparse_namespace(**kwargs):
    from argparse import Namespace
    return Namespace(**kwargs)


if __name__=='__main__':unittest.main()
