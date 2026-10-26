"""Fixture-only orchestration checks: no model execution or live qualification.

The full-route fixtures mock native preflight/view/model endpoints. Only pinned
upstream prepare/extract and ordinary bgzip/tabix/chromosome renaming run on tiny
temporary inputs. Fixture profiles cannot qualify or enable production.
"""
import argparse
import contextlib
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
SPEC = importlib.util.spec_from_file_location('cadd_shared_workflow_fixture', ROOT / 'scripts/cadd_shared_workflow.py')
M = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(M)
Y = Path('/stornext/snfs170/next-gen/scratch/yangyxt')
D = Path('/stornext/snfs170/next-gen/scratch/yangyxt_dt/projects/PriVA_upgrade')
UPSTREAM = Y / 'git_repos/CADD-scripts-v1.7.2-review'
PYTHON2 = D / 'caches/runtime/cadd-1.7.2/bin/cadd-python2'
BIN = Y / 'mamba/envs/priva_acmg/bin'


class SharedWorkflowTests(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory(prefix='cadd-workflow-fixture-', dir=D / 'tmp')
        self.root = Path(self.temporary.name)
        self.scratch = self.root / 'scratch'; self.scratch.mkdir()
        self.environment = mock.patch.dict(os.environ, PATH=str(BIN) + ':' + os.environ['PATH'],
                                           TMPDIR=str(self.scratch), PYTHONDONTWRITEBYTECODE='1')
        self.environment.start()

    def tearDown(self):
        self.environment.stop()
        self.temporary.cleanup()

    def row(self, pos='10', alt='G', payload='same', chrom='1', ref='A'):
        return [chrom, pos, '.', ref, alt, '.', '.', 'AF=0.1,0.2;EsmScoreMissense=stale;CSQ=' + payload, 'GT', '0/1']

    def write_vcf(self, path, rows):
        text = ('##fileformat=VCFv4.2\n##contig=<ID=chr1,length=1000>\n'
                '##INFO=<ID=AF,Number=A,Type=Float,Description="frequency">\n'
                '##INFO=<ID=EsmScoreMissense,Number=.,Type=String,Description="stale fixture model tag">\n'
                '##INFO=<ID=CSQ,Number=.,Type=String,Description="Format: Allele|Feature">\n'
                '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n'
                '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tsample\n')
        text += ''.join('\t'.join(row) + '\n' for row in rows)
        with path.open('wb') as output:
            subprocess.run([str(BIN / 'bgzip'), '-c'], input=text.encode(), stdout=output, check=True, timeout=10)
        return path

    def disk_mapping(self, original, prepared, converted=None):
        folder=Path(tempfile.mkdtemp(prefix='mapping-',dir=self.scratch))
        headers='##fileformat=VCFv4.2\n##INFO=<ID=CSQ,Number=.,Type=String,Description="Format: Allele|Feature">\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n'
        source=folder/'source.vcf'; renamed=folder/'renamed.vcf'; authority=folder/'prepared.vcf'; proof=folder/'mapping.jsonl'
        source.write_text(headers+''.join('\t'.join(row)+'\n' for row in original))
        if converted is None:
            converted=[['1' if row[0]=='chr1' else row[0], *row[1:]] for row in original]
        renamed.write_text(headers+''.join('\t'.join(row)+'\n' for row in converted))
        authority.write_text(headers+''.join('\t'.join(row)+'\n' for row in prepared))
        return M.OccurrenceMap(source,renamed,authority,{'chr1':'1'},folder,proof), authority, proof

    def test_reorder_and_equivalent_duplicate_collapse_keep_all_occurrences(self):
        rows = [self.row('20', 'T', 'twenty'), self.row(), self.row()]
        prepared = [['1', '10', '.', 'A', 'G'], ['1', '20', '.', 'A', 'T']]
        mapper, authority, proof_path = self.disk_mapping(rows,prepared)
        with mapper: mapped=list(mapper.shared_records(authority))
        proof=[json.loads(line) for line in proof_path.read_text().splitlines() if json.loads(line)['stage']=='prepared_occurrence']
        self.assertEqual([M.key(row) for row in mapped], [M.key(row) for row in prepared])
        self.assertEqual([item['source_occurrence'] for item in proof if item['prepared_occurrence']==0], [1, 2])
        self.assertEqual([item['source_occurrence'] for item in proof if item['prepared_occurrence']==1], [0])
        self.assertEqual(mapped[0], ['1', '10', '.', 'A', 'G', '.', '.', 'CSQ=same'])

    def test_ambiguous_duplicate_and_changed_prepared_representation_fail(self):
        with self.assertRaisesRegex(M.WorkflowError, 'differing complete CSQs'):
            self.disk_mapping([self.row(payload='one'), self.row(payload='two')], [['1', '10', '.', 'A', 'G']])
        with self.assertRaisesRegex(M.WorkflowError, 'changed/dropped allele representation'):
            self.disk_mapping([self.row('100', 'AAT', ref='AAC')], [['1', '102', '.', 'C', 'T']])

    def test_private_novel_vcf_has_only_first_csq_info_and_no_stale_metadata(self):
        headers = ['##fileformat=VCFv4.2', '##INFO=<ID=AF,Number=A,Type=Float,Description="AF">',
                   '##INFO=<ID=CSQ,Number=.,Type=String,Description="Format: Allele|Feature">', '#CHROM\tPOS']
        path = self.root / 'private.gz'
        mapper, authority, _ = self.disk_mapping([self.row(payload='G|T1,G|T2')], [['1', '10', '.', 'A', 'G']])
        with mapper: M.write_shared(path, headers, mapper.shared_records(authority))
        output_headers, output_rows = M.vcf(path)
        output_rows=list(output_rows)
        infos = [line for line in output_headers if line.startswith('##INFO=')]
        self.assertEqual(len(infos), 1); self.assertTrue(infos[0].startswith('##INFO=<ID=CSQ,'))
        self.assertEqual(len(output_rows[0]), 8)
        self.assertEqual(output_rows[0][7], 'CSQ=G|T1,G|T2')
        self.assertNotIn('EsmScore', gzip.open(path, 'rt').read())

    def test_chromosome_rename_proves_all_fields_and_occurrence_permutation(self):
        original = [self.row('20', chrom='chr1'), self.row(chrom='chr1')]
        converted = [self.row(), self.row('20')]
        mapper, _, proof=self.disk_mapping(original,[['1','10','.','A','G'],['1','20','.','A','G']],converted)
        mapper.close()
        self.assertEqual([json.loads(line)['source_occurrence'] for line in proof.read_text().splitlines()
                          if json.loads(line)['stage']=='chromosome_rename'],[1,0])
        converted[0][7] = 'CSQ=changed'
        with self.assertRaisesRegex(M.WorkflowError, 'changed records'):
            self.disk_mapping(original,[['1','10','.','A','G'],['1','20','.','A','G']],converted)

    def test_publisher_split_keeps_every_matching_feature_row_but_no_other_alt(self):
        previous = [['1','10','.','A','G'], ['1','10','.','A','T']]
        M.prove_split(previous, [previous[1]], ['1\t10\tA\tG\tT1', '1\t10\tA\tG\tT2'])
        with self.assertRaisesRegex(M.WorkflowError, 'unrelated allele'):
            M.prove_split(previous, [previous[1]], ['1\t10\tA\tC\tbad'])
        with self.assertRaisesRegex(M.WorkflowError, 'split membership or order'):
            M.prove_split(previous, [], ['1\t10\tA\tG\tT1'])

    def test_original_join_last_zero_hit_header_and_whole_line_ties_keep_duplicates(self):
        pre, novel = self.root/'pre.tsv', self.root/'novel.tsv'
        pre.write_text('## metadata\n#first\n1\t10\tA\tG\tT2\n1\t10\tA\tG\tT1\n1\t10\tA\tG\tT1\n#last-zero-hit\n')
        output = self.root/'joined.tsv'
        M.join_scores([pre], output, dict(os.environ, LC_ALL='C'))
        self.assertEqual(output.read_text(), '#last-zero-hit\n1\t10\tA\tG\tT1\n1\t10\tA\tG\tT1\n1\t10\tA\tG\tT2\n')
        novel.write_text('#novel-header\n1\t2\tA\tG\tnew\n')
        output2=self.root/'mixed.tsv'; M.join_scores([pre, novel], output2, dict(os.environ, LC_ALL='C'))
        self.assertTrue(output2.read_text().startswith('#novel-header\n1\t2\tA\tG\tnew\n'))
        empty=self.root/'empty.tsv'; empty.write_text('## no actual header\n')
        with self.assertRaisesRegex(M.WorkflowError, 'no score header'):
            M.join_scores([empty], self.root/'must-not-exist', dict(os.environ))

    def cache(self, name, lines):
        path = self.root / (name + '.tsv.gz')
        header='#Chrom\tPos\tRef\tAlt\tFeatureID\tRawScore\tPHRED\n'
        with path.open('wb') as output:
            subprocess.run([str(BIN/'bgzip'), '-c'], input=(header+''.join(line+'\n' for line in lines)).encode(), stdout=output, check=True)
        subprocess.run([str(BIN/'tabix'), '-s','1','-b','2','-e','2',str(path)], check=True)
        return path

    def execution_fixture(self, rows, caches, *, missing_last_index=False):
        input_path=self.write_vcf(self.root/'input.vcf.gz', rows)
        config={'tmp_dir':str(self.scratch), 'assembly':'GRCh38'}
        profile_path=self.root/'fixture-profile.json'
        profile_path.write_text('{"fixture_only":true,"not_live_qualification":true}\n')
        profile={'contract':{'columns':['Chrom','Pos','Ref','Alt','FeatureID'], 'resource_manifest_sha256':'f'*64},
                 'annotation':{'env':{'CADD_PYTHON2':str(PYTHON2)}}}
        block={'prepare':{'path':str(UPSTREAM/'src/scripts/VCF2vepVCF.py')},
               'extract':{'path':str(UPSTREAM/'src/scripts/extract_scored.py')},
               'files':[{'path':str(cache)} for cache in caches]}
        task, shared, view, native = (mock.Mock() for _ in range(4))
        captured=[]
        def build(source, profile, output):
            headers, records=M.vcf(source); captured.append((headers,list(records)))
            output.mkdir(); M.save(output/'manifest.json', {'fixture_only':True})
        view.build.side_effect=build
        def run(config, profile_path, bundle, output):
            output.mkdir(); calculation=output/'calculation_bundle'; calculation.mkdir()
            _, novel_rows=captured[-1]
            (calculation/'scores.tsv').write_text('#Chrom\tPos\tRef\tAlt\tFeatureID\tRawScore\tPHRED\n'+''.join(
                '\t'.join([row[0],row[1],row[3],row[4],'NOVEL','0.2','7'])+'\n' for row in novel_rows))
        native.run.side_effect=run
        if missing_last_index: Path(str(caches[-1])+'.tbi').unlink()
        args=argparse.Namespace(config=str(self.root/'fixture.yaml'), profile=str(profile_path), input=str(input_path), output=str(self.root/'output.tsv'))
        args_config=(config, profile_path, profile, block, {}, (task, shared, view, native))
        return args, args_config, native, captured

    def execute(self, args, setup):
        with mock.patch.object(M,'preflight',return_value=setup), mock.patch.object(M,'native_rows'), \
             mock.patch.object(M,'adapter_identity',return_value={'fixture_only':'no production source'}):
            return M.calculate(args)

    def test_actual_original_all_prescored_multiple_rows_first_file_wins_no_model(self):
        a=self.cache('A', ['1\t10\tA\tG\tT2\t0.1\t5','1\t10\tA\tG\tT1\t0.1\t5'])
        b=self.cache('B', ['1\t10\tA\tG\tLATER\t99\t99'])
        args, setup, native, captured=self.execution_fixture([self.row(chrom='chr1')], [a,b])
        self.execute(args,setup)
        native.run.assert_not_called(); self.assertEqual(captured, [])
        _, rows=M.table(args.output)
        rows=list(rows)
        self.assertEqual(len(rows),2); self.assertTrue(all('LATER' not in row for row in rows))
        self.assertTrue(rows[0].startswith('1\t10\tA\tG\tT1'))
        receipt=M.read_json(args.output+'.cadd-native.json')
        self.assertEqual(receipt['publisher_files_processed'],2); self.assertEqual(receipt['novel_alleles'],0)
        self.execute(args,setup); native.run.assert_not_called()

    def test_actual_original_late_missing_index_after_all_hits_still_fails(self):
        a=self.cache('A',['1\t10\tA\tG\tHIT\t0.1\t5']); b=self.cache('B',[])
        args,setup,native,_=self.execution_fixture([self.row(chrom='chr1')],[a,b],missing_last_index=True)
        with self.assertRaises(subprocess.CalledProcessError) as failure: self.execute(args,setup)
        self.assertIn('extract_scored.py', failure.exception.cmd[1])
        native.run.assert_not_called(); self.assertFalse(Path(args.output).exists())

    def test_plain_and_gzipped_sources_keep_bytes_and_produce_identical_plain_scores(self):
        cache=self.cache('A',['1\t10\tA\tG\tT2\t0.1\t5','1\t10\tA\tG\tT1\t0.1\t5'])
        args,setup,native,_=self.execution_fixture([self.row(chrom='chr1')],[cache])
        compressed_before=Path(args.input).read_bytes()
        self.execute(args,setup); expected=Path(args.output).read_bytes()
        self.assertEqual(Path(args.input).read_bytes(),compressed_before)
        self.assertTrue(Path(args.output+'.native/source.original.vcf.gz').is_file())
        plain=self.root/'plain.vcf'
        plain.write_text(gzip.open(args.input,'rt').read()); plain_before=plain.read_bytes()
        plain_args=argparse.Namespace(**vars(args)); plain_args.input=str(plain); plain_args.output=str(self.root/'plain-output.tsv')
        self.execute(plain_args,setup)
        self.assertEqual(Path(plain_args.output).read_bytes(),expected)
        self.assertEqual(plain.read_bytes(),plain_before)
        self.assertEqual(Path(plain_args.output+'.native/source.original.vcf').read_bytes(),plain_before)
        native.run.assert_not_called()

    def test_mixed_hits_and_novel_uses_only_exact_allele_with_clean_shared_info(self):
        a=self.cache('A',['1\t10\tA\tG\tHIT\t0.1\t5'])
        args,setup,native,captured=self.execution_fixture([self.row('20','T',chrom='chr1'), self.row(chrom='chr1')],[a])
        original=Path(args.input).read_bytes(); self.execute(args,setup)
        native.run.assert_called_once(); self.assertEqual(Path(args.input).read_bytes(),original)
        headers, rows=captured[0]
        self.assertEqual([M.key(row) for row in rows],[('1','20','A','T')])
        self.assertEqual(rows[0][7],'CSQ=same'); self.assertEqual(len(rows[0]),8)
        self.assertEqual(len([header for header in headers if header.startswith('##INFO=')]),1)
        _, final=M.table(args.output); self.assertEqual(len(list(final)),2)

    def test_unreceipted_warm_output_is_preserved_and_not_accepted(self):
        a=self.cache('A',[]); args,setup,native,_=self.execution_fixture([self.row(chrom='chr1')],[a])
        Path(args.output).write_text('historical untouched output\n')
        with self.assertRaisesRegex(M.WorkflowError,'no matching closed native receipt'): self.execute(args,setup)
        self.assertEqual(Path(args.output).read_text(),'historical untouched output\n'); native.run.assert_not_called()
        self.assertFalse(Path(args.output+'.native').exists())

    def test_unqualified_preflight_fails_before_any_output_or_model_work(self):
        args=argparse.Namespace(config='fixture',profile=None,input='not-read',output=str(self.root/'new.tsv'))
        with mock.patch.object(M,'preflight',side_effect=M.WorkflowError('unqualified')), \
             mock.patch.object(M,'digest') as digest, self.assertRaisesRegex(M.WorkflowError,'unqualified'):
            M.calculate(args)
        digest.assert_not_called(); self.assertEqual(list(self.root.iterdir()),[self.scratch])

    def test_missing_explicit_profile_never_selects_installed_candidates(self):
        candidates=self.root/'runtime/cadd-1.7.2/profiles'; candidates.mkdir(parents=True)
        for name in ('GRCh38.installed.json','GRCh38.vep116.fixture.installed.json'):
            (candidates/name).write_text('{"status":"installed_not_scientifically_qualified"}\n')
        before={path.name:path.read_bytes() for path in candidates.iterdir()}
        with self.assertRaisesRegex(M.WorkflowError,'--cadd-profile or config cadd_profile') as error:
            M.selection(argparse.Namespace(profile=None),{'cache_root':str(self.root),'assembly':'GRCh38'})
        self.assertIn('CADD installation receipt',str(error.exception))
        self.assertNotIn(str(candidates),str(error.exception))
        self.assertEqual({path.name:path.read_bytes() for path in candidates.iterdir()},before)

    def test_same_position_alt_miss_extraction_preserves_original_exact_rows(self):
        input_path=self.write_vcf(self.root/'source.gz',[self.row(alt='G',chrom='chr1'),self.row(alt='T',chrom='chr1')])
        covered=self.root/'covered.tsv'; covered.write_text('#Chrom\tPos\tRef\tAlt\tFeatureID\tPHRED\nchr1\t10\tA\tG\tHIT\t5\n')
        config=self.root/'fixture.yaml'; config.write_text('tmp_dir: '+str(self.scratch)+'\n')
        output=self.root/'miss.vcf.gz'
        args=argparse.Namespace(config=str(config),profile=None,input=str(input_path),covered=str(covered),output=str(output))
        with mock.patch.object(M,'preflight'): M.extract_misses(args)
        _, rows=M.vcf(output)
        rows=list(rows)
        self.assertEqual([M.key(row) for row in rows],[('chr1','10','A','T')])
        self.assertEqual(rows[0][7],self.row(alt='T',chrom='chr1')[7])

    def shell_fixture(self, purpose, *, existing_csq=False):
        scripts=self.root/'shell-scripts'; scripts.mkdir()
        binary=self.root/'shell-bin'; binary.mkdir()
        (scripts/'argparse.bash').write_text('argparse() { while [[ $# -gt 0 ]]; do case "$1" in --input_vcf) input_vcf=$2;; --config) config=$2;; --threads) threads=$2;; --cadd-profile) cadd_profile=$2;; esac; shift 2; done; }\n')
        for name in ('anno_vep_args',): (scripts/name).write_text('# fixture parser uses named inputs\n')
        helper=scripts/'cadd_shared_workflow.py'
        helper.write_text('import sys\nif sys.argv[1]=="vep-options": print("/fixture/plugins\\nEnsembl\\n116\\n0")\nelif sys.argv[1]=="verify-source" and '+repr(existing_csq)+': raise SystemExit(2)\n')
        (scripts/'clinvar_gene_refresh.py').write_text('# fixture-only installer receipt endpoint\n')
        (self.root/'input.vcf').write_text('##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n1\t10\t.\tA\tG\t.\t.\t.\n')
        (self.root/'input.vcf.gz').write_bytes(b'fixture; subprocess VCF tools are mocks')
        (self.root/'fixture.yaml').write_text('assembly: GRCh38\n')
        program=binary/'vep'
        program.write_text('#!'+str(BIN/'python')+'\nimport json,os,pathlib,shutil,sys\na=sys.argv[1:]\npathlib.Path(os.environ["FIXTURE"],"vep-argv.json").write_text(json.dumps({"argv":a,"hash_seed":os.environ.get("PERL_HASH_SEED")}))\nshutil.copyfile(pathlib.Path(os.environ["FIXTURE"],"input.vcf"),a[a.index("-o")+1])\n')
        program.chmod(0o755)
        script='''source "$ROOT/scripts/annotation_vcf.sh"
SCRIPT_DIR="$FIXTURE/shell-scripts"
read_yaml() { case "$2" in assembly) printf GRCh38;; threads) printf 2;; pext_enabled) printf false;; *) printf '%s/file' "$FIXTURE";; esac; }
priva_require_legacy_spliceai() { return 0; }
check_path() { return 0; }
check_vcf_infotags() { return "$EXISTING_CSQ_STATUS"; }
display_vcf() { return 0; }
announce_remove_tmps() { return 0; }
bcftools() { local previous= output= item; for item in "$@"; do [[ $previous != -o ]] || output=$item; previous=$item; done; cp -- "${@: -1}" "$output"; }
tabix() { touch -- "${@: -1}.tbi"; }
"$PURPOSE" --input_vcf "$FIXTURE/input.vcf.gz" --config "$FIXTURE/fixture.yaml" --threads 2 "${CADD_OPTIONS[@]}"
'''
        # Bash arrays must be constructed by trusted fixture code, never eval'd
        # from a profile/configuration string.
        prefix='CADD_OPTIONS=(--cadd-profile /fixture/qualified.json)\n' if purpose=='anno_VEP_data' else 'CADD_OPTIONS=()\n'
        env=dict(os.environ,ROOT=str(ROOT),FIXTURE=str(self.root),PURPOSE=purpose,
                 EXISTING_CSQ_STATUS='0' if existing_csq else '1',PATH=str(binary)+':'+os.environ['PATH'])
        return subprocess.run(['bash','-c',prefix+script],env=env,capture_output=True,text=True,timeout=15)

    def test_runtime_vep_adds_qualified_selection_without_filtering_priva_csqs(self):
        result=self.shell_fixture('anno_VEP_data')
        self.assertEqual(result.returncode,0,result.stderr)
        observed=M.read_json(self.root/'vep-argv.json'); argv=observed['argv']
        for flag in ('--merged','--use_transcript_ref','--regulatory','--ccds','--domains','--numbers','--canonical','--total_length'):
            self.assertIn(flag,argv)
        self.assertEqual(argv[argv.index('--sift')+1],'b'); self.assertEqual(argv[argv.index('--polyphen')+1],'b')
        self.assertIn('CaddSelection,cadd_transcript_sources=Ensembl',argv)
        self.assertTrue(any(item.startswith('AlphaMissense,') for item in argv))
        self.assertTrue(any(item.startswith('SpliceAI,') for item in argv))
        self.assertNotIn('--per_gene',argv); self.assertNotIn('--pick',argv)
        self.assertEqual(observed['hash_seed'],'0')

    def test_installer_candidate_does_not_activate_cadd_options(self):
        result=self.shell_fixture('anno_VEP_cache_candidate')
        self.assertEqual(result.returncode,0,result.stderr)
        argv=M.read_json(self.root/'vep-argv.json')['argv']
        self.assertNotIn('--regulatory',argv); self.assertNotIn('--sift',argv)
        self.assertFalse(any('CaddSelection' in item for item in argv))

    def test_old_csq_input_is_rejected_without_second_vep(self):
        result=self.shell_fixture('anno_VEP_data',existing_csq=True)
        self.assertNotEqual(result.returncode,0)
        self.assertFalse((self.root/'vep-argv.json').exists())

    def test_runtime_options_log_pinned_profile_and_native_context(self):
        profile_path=self.root/'logging-fixture.json'; profile_path.write_text('{"fixture_only":true}\n')
        shared=mock.Mock(); shared.verify_profile.return_value=({'selection_context':{
            'cadd_transcript_sources':['Ensembl'], 'perl_hash_seed':'0', 'annotation_options':{'regulatory':'1'}},
            'vep_release':'116'}, 'a'*64)
        fixture=({},profile_path,{},None,None,(None,shared,None,None))
        output,error=io.StringIO(),io.StringIO()
        with mock.patch.object(M,'preflight',return_value=fixture), contextlib.redirect_stdout(output), contextlib.redirect_stderr(error):
            M.vep_options(argparse.Namespace())
        self.assertEqual(len(output.getvalue().splitlines()),4)
        logged=json.loads(error.getvalue())
        self.assertEqual(logged['profile_sha256'],M.digest(profile_path))
        self.assertEqual(logged['selection_context_sha256'],'a'*64)
        self.assertEqual(logged['annotation_options'],{'regulatory':'1'})


if __name__ == '__main__': unittest.main()
