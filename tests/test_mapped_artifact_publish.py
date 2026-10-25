"""Tiny generated pickles and real wrapper calls; no scientific calculation."""
import importlib.util
import json
import os
from pathlib import Path
import pickle
import shutil
import subprocess
import sys
import tempfile
from types import SimpleNamespace
import unittest
from unittest import mock


ROOT = Path(__file__).resolve().parents[1]
SPEC = importlib.util.spec_from_file_location("mapped_publish_test", ROOT / "scripts/mapped_artifact_publish.py")
M = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(M)
PAYLOAD = {"ordered_rows": [("B", 0.4), ("A", None), ("B", 0.4)], "regions": {"ENST1": [(3, 5)]}}
PRODUCER = '''import json,os,pickle,signal,sys
from pathlib import Path
a=sys.argv[1:]
with open(os.environ['ARGV_LOG'],'a') as f: f.write(json.dumps(a)+'\\n')
p=Path(a[a.index('--output')+1])
if os.environ.get('FAIL_PRODUCER'):
    p.write_bytes(b'incomplete generated pickle')
    if os.environ['FAIL_PRODUCER']=='signal': os.kill(os.getpid(),signal.SIGTERM)
    sys.exit(9)
p.write_bytes(pickle.dumps({'ordered_rows':[('B',0.4),('A',None),('B',0.4)],'regions':{'ENST1':[(3,5)]}},protocol=pickle.HIGHEST_PROTOCOL))
'''


class PublicationFixtures(unittest.TestCase):
    def setUp(self):
        tmp = tempfile.TemporaryDirectory(prefix="mapped_publish_", dir=os.environ.get("PRIVA_TEST_TMPDIR", os.environ.get("TMPDIR")))
        self.addCleanup(tmp.cleanup)
        self.folder = Path(tmp.name)
        self.cache = self.folder / "managed cache"
        self.output_dir = self.cache / "plugins/AlphaMissense"
        self.output_dir.mkdir(parents=True)
        self.scripts = self.folder / "scripts"
        self.scripts.mkdir()
        for name in ("mapped_artifact_publish.py", "biomart_nm_to_enst.py", "biomart_ftp_nm_to_enst.py"):
            shutil.copy2(ROOT / "scripts" / name, self.scripts / name)
        for name in M.PRODUCERS.values():
            (self.scripts / name).write_text(PRODUCER)
        self.mapping = self.folder / "selected mapping.tsv"
        self.mapping.write_text("selected mapping fixture\n")
        self.anchors = []
        for relative in ("gnomAD_RMC/gnomad_v2.1.1_rmc_flat.hg38.tsv", "HCSeeker/HC_spots.hg38.tsv"):
            path = self.folder / "data" / relative
            path.parent.mkdir(parents=True, exist_ok=True)
            path.write_text("fixture anchor\n")
            self.anchors.append(path)
        self.config = self.folder / "runtime config.json"
        self.config.write_text('{"fixture":"actual runtime config"}\n')
        self.argv_log = self.folder / "argv.jsonl"
        env = mock.patch.dict(os.environ, {"ARGV_LOG": str(self.argv_log), "PYTHONDONTWRITEBYTECODE": "1"})
        env.start()
        self.addCleanup(env.stop)

    def args(self, kind="motifs", target=None):
        target = target or self.output_dir / ("motifs.pkl" if kind == "motifs" else "pm1.pkl")
        return SimpleNamespace(artifact=target, mapping=self.mapping, input=self.anchors,
                               cache_root=self.cache, output_dir=self.output_dir, assembly="hg38",
                               config=self.config, producer=self.scripts / M.PRODUCERS[kind], kind=kind,
                               expected_identity=json.dumps(M.BIOMART.mapping_identity(self.mapping, self.anchors)),
                               arguments=["--literal", "space and ' quote", "--output", str(target), "--skip_das"])

    def complete(self, args, content=b"previous completed pickle"):
        args.artifact.write_bytes(content)
        M.BIOMART.atomic_json(M.receipt_path(args.artifact),
                              {"mapping_identity": M.BIOMART.mapping_identity(args.mapping, args.input),
                               "artifact_bytes": len(content), "artifact_sha256": M.BIOMART.sha256(args.artifact)})

    def pair(self, target):
        return [M.file_state(path) for path in (target, M.receipt_path(target))]

    def interrupt(self, args, phase):
        real_replace = os.replace
        fired = False

        def replace(source, destination):
            nonlocal fired
            result = real_replace(source, destination)
            dest = Path(destination)
            selected = {"marker": M.pending_path(args.artifact), "artifact": args.artifact,
                        "receipt": M.receipt_path(args.artifact)}
            hit = dest == selected.get(phase)
            if phase == "staged_receipt":
                hit = dest.name.endswith(".nm_to_enst.json") and dest.parent != args.artifact.parent
            if hit and not fired:
                fired = True
                raise OSError("injected interruption after " + phase)
            return result

        with mock.patch.object(M.os, "replace", side_effect=replace):
            with self.assertRaisesRegex(OSError, "injected interruption"):
                M.run(args)
        self.assertTrue(fired)


class PublicationTransactions(PublicationFixtures):
    def test_run_changes_only_output_and_preserves_exact_pickle_and_receipt(self):
        for kind in M.PRODUCERS:
            with self.subTest(kind=kind):
                args = self.args(kind)
                before = list(args.arguments)
                M.run(args)
                observed = json.loads(self.argv_log.read_text().splitlines()[-1])
                output_index = before.index("--output") + 1
                self.assertNotEqual(observed[output_index], before[output_index])
                observed[output_index] = before[output_index]
                self.assertEqual(observed, before)
                self.assertEqual(args.artifact.read_bytes(), pickle.dumps(PAYLOAD, protocol=pickle.HIGHEST_PROTOCOL))
                M.BIOMART.check_receipt(args.artifact, self.mapping, self.anchors)
                self.assertFalse(M.pending_path(args.artifact).exists())

    def test_failed_or_signalled_pickle_write_preserves_existing_pair(self):
        for failure in ("exit", "signal"):
            with self.subTest(failure=failure):
                args = self.args(target=self.output_dir / (failure + ".pkl"))
                self.complete(args)
                before = self.pair(args.artifact)
                with mock.patch.dict(os.environ, {"FAIL_PRODUCER": failure}):
                    with self.assertRaises(subprocess.CalledProcessError):
                        M.run(args)
                self.assertEqual(self.pair(args.artifact), before)
                self.assertFalse(M.pending_path(args.artifact).exists())
                self.assertTrue(list(self.output_dir.glob(".mapped-artifact-motifs-*/" + args.artifact.name)))

    def test_interruption_during_staged_receipt_never_exposes_new_pickle(self):
        args = self.args()
        self.complete(args)
        before = self.pair(args.artifact)
        self.interrupt(args, "staged_receipt")
        self.assertEqual(self.pair(args.artifact), before)
        self.assertFalse(M.pending_path(args.artifact).exists())
        M.run(args)  # Failed private stages do not block a fresh calculation.
        M.BIOMART.check_receipt(args.artifact, self.mapping, self.anchors)

    def test_each_commit_interruption_recovers_without_rerunning_producer(self):
        for kind in M.PRODUCERS:
            for phase in ("marker", "artifact", "receipt"):
                with self.subTest(kind=kind, phase=phase):
                    args = self.args(kind, self.output_dir / (kind + "-" + phase + ".pkl"))
                    self.complete(args)
                    self.interrupt(args, phase)
                    calls = self.argv_log.read_bytes()
                    self.assertTrue(M.recover(args))
                    self.assertEqual(self.argv_log.read_bytes(), calls)
                    self.assertFalse(M.recover(args))
                    self.assertEqual(pickle.loads(args.artifact.read_bytes()), PAYLOAD)
                    M.BIOMART.check_receipt(args.artifact, self.mapping, self.anchors)

    def test_changed_mapping_pending_namespace_recovers_before_normal_routing(self):
        for kind in M.PRODUCERS:
            with self.subTest(kind=kind):
                original = self.args(kind)
                self.complete(original)
                before = self.pair(original.artifact)
                self.mapping.write_text("changed mapping " + kind + "\n")
                target = M.BIOMART.artifact_route(artifact=original.artifact, mapping=self.mapping,
                    cache_root=self.cache, output_dir=self.output_dir, assembly="hg38", kind=kind, inputs=self.anchors)
                generated = self.args(kind, target)
                self.interrupt(generated, "artifact")
                self.assertTrue(M.recover(self.args(kind, original.artifact)))
                self.assertEqual(self.pair(original.artifact), before)
                self.assertEqual(M.BIOMART.artifact_route(artifact=original.artifact, mapping=self.mapping,
                    cache_root=self.cache, output_dir=self.output_dir, assembly="hg38", kind=kind, inputs=self.anchors), target)

    def test_unreceipted_external_or_orphan_receipt_is_never_recovery_authority(self):
        for orphan in (False, True):
            with self.subTest(orphan=orphan):
                args = self.args(target=self.output_dir / (str(orphan) + ".pkl"))
                path = M.receipt_path(args.artifact) if orphan else args.artifact
                path.write_bytes(b"unknown supplied bytes")
                before = M.file_state(path)
                self.assertFalse(M.recover(args))
                with self.assertRaises(ValueError):
                    M.run(args)
                self.assertEqual(M.file_state(path), before)
                self.assertFalse(self.argv_log.exists())

    def test_changed_producer_config_mapping_prior_or_staged_content_refuses_recovery(self):
        for mutation in ("producer", "config", "mapping", "prior", "staged"):
            with self.subTest(mutation=mutation):
                args = self.args(target=self.output_dir / (mutation + ".pkl"))
                self.complete(args)
                self.interrupt(args, "marker")
                marker = json.loads(M.pending_path(args.artifact).read_text())
                changed = {"producer": args.producer, "config": args.config,
                           "mapping": args.mapping, "prior": args.artifact,
                           "staged": Path(marker["stage"]) / args.artifact.name}[mutation]
                content = changed.read_bytes()
                changed.write_bytes(content + b"changed")
                before = self.pair(args.artifact)
                with self.assertRaises(ValueError):
                    M.recover(args)
                self.assertEqual(self.pair(args.artifact), before)
                self.assertTrue(M.pending_path(args.artifact).exists())
                changed.write_bytes(content)

    def test_malformed_linked_and_escaped_markers_refuse_without_mutation(self):
        for mutation in ("format", "stage_escape", "target_escape", "marker_link"):
            with self.subTest(mutation=mutation):
                args = self.args(target=self.output_dir / (mutation + ".pkl"))
                self.complete(args)
                self.interrupt(args, "marker")
                marker = M.pending_path(args.artifact)
                value = json.loads(marker.read_text())
                if mutation == "format":
                    marker.write_text('{}')
                elif mutation == "stage_escape":
                    value["stage"] = str(self.folder)
                    marker.write_text(json.dumps(value))
                elif mutation == "target_escape":
                    value["target"] = str(self.folder / "outside.pkl")
                    marker.write_text(json.dumps(value))
                else:
                    saved = self.folder / "marker-target.json"
                    marker.rename(saved)
                    marker.symlink_to(saved)
                before = self.pair(args.artifact)
                with self.assertRaises(ValueError):
                    M.recover(args)
                self.assertEqual(self.pair(args.artifact), before)

    def test_receipt_first_impossible_state_is_rejected(self):
        args = self.args()
        self.complete(args)
        self.interrupt(args, "marker")
        value = json.loads(M.pending_path(args.artifact).read_text())
        os.replace(M.receipt_path(Path(value["stage"]) / args.artifact.name), M.receipt_path(args.artifact))
        before = self.pair(args.artifact)
        with self.assertRaisesRegex(ValueError, "ordered two-file"):
            M.recover(args)
        self.assertEqual(self.pair(args.artifact), before)

    def test_private_stage_allows_inherited_setgid_but_rejects_group_access(self):
        for mode in (0o2700, 0o750):
            with self.subTest(mode=oct(mode)):
                args = self.args(target=self.output_dir / (str(mode) + ".pkl"))
                self.complete(args)
                self.interrupt(args, "marker")
                value = json.loads(M.pending_path(args.artifact).read_text())
                Path(value["stage"]).chmod(mode)
                before = self.pair(args.artifact)
                if mode == 0o2700:
                    self.assertTrue(M.recover(args))
                else:
                    with self.assertRaisesRegex(ValueError, "stage is not owned"):
                        M.recover(args)
                    self.assertEqual(self.pair(args.artifact), before)

    def test_linked_stage_payload_is_rejected_before_receipt_creation(self):
        args = self.args()
        outside = self.folder / "outside.pkl"
        outside.write_bytes(b"supplied bytes")
        before = M.file_state(outside)
        def child(argv, **kwargs):
            Path(argv[argv.index("--output") + 1]).symlink_to(outside)
        with mock.patch.object(M.subprocess, "run", side_effect=child):
            with self.assertRaisesRegex(ValueError, "regular"):
                M.run(args)
        self.assertFalse(args.artifact.exists())
        self.assertEqual(M.file_state(outside), before)


class InstallerWrappers(PublicationFixtures):
    def setup_wrapper(self, kind):
        # This labelled stub stands only for unrelated mapping prerequisites.
        (self.scripts / "installer_prerequisites.py").write_text(
            "import json,sys\nprint(json.load(open(sys.argv[1]))['nm_to_enst_tsv'])\n")
        # This suite isolates mapping/pair publication. Complete-base science
        # and native migration have real BGZF/receipt coverage in the separate
        # test_alphamissense_annotation_base suite.
        (self.scripts / "alphamissense_annotation_base.py").write_text('''import json,sys
from pathlib import Path
def motif_binding(*args): return {"fixture":"already verified complete base"}
if __name__ == "__main__":
    a=sys.argv[1:]
    if a[0] == "path": print(json.load(open(a[a.index("--config")+1]))["alphamissense_vep_vcf"])
    elif a[0] == "motif-check":
        r=Path(a[a.index("--artifact")+1]+".nm_to_enst.json")
        sys.exit(0 if r.exists() and json.loads(r.read_text()).get("scientific_input")==motif_binding() else 1)
''')
        artifact = self.output_dir / ("motifs.pkl" if kind == "motifs" else "pm1.pkl")
        source_vcf = self.cache / "input.vcf.gz"
        source_vcf.write_bytes(b"not read by this publication fixture")
        domains = self.cache / "all_intolerant.pkl"
        domains.write_bytes(b"not read by this publication fixture")
        self.config.write_text(json.dumps({"assembly": "GRCh38", "threads": 2, "cache_root": str(self.cache),
            "vep_plugins_cachedir": str(self.cache / "plugins"), "nm_to_enst_tsv": str(self.mapping),
            "alphamissense_intolerant_motifs": str(self.output_dir / "motifs.pkl"),
            "pm1_regions_pkl": str(self.output_dir / "pm1.pkl"), "all_intolerant_domains": str(domains),
            "alphamissense_pd_stat": str(self.cache / "domain.stats.pkl"),
            "alphamissense_vep_vcf": str(source_vcf)}))
        args = self.args(kind, artifact)
        if kind == "motifs":
            args.annotation_base = source_vcf
        return args

    def wrapper(self, kind):
        function = "AlphaMissense_pick_missense_constrained_segments" if kind == "motifs" else "prepare_pm1_regions"
        body = r'''
source "$1"
SCRIPT_DIR=$2
BASE_DIR=$3
DATA_DIR="$3/data"
log() { printf '%s\n' "$*" >&2; }
read_yaml() { python -c 'import json,sys;v=json.load(open(sys.argv[1])).get(sys.argv[2]);print("null" if v is None else v)' "$@"; }
update_yaml() { [[ ${FAIL_CONFIG_WRITE:-0} != 1 ]] || return 17; python -c 'import json,sys;p,k,v=sys.argv[1:];d=json.load(open(p));d[k]=v;open(p,"w").write(json.dumps(d))' "$@"; }
yaml_value_or_default() { local v;v=$(read_yaml "$1" "$2"); if [[ $v == null || -z $v ]];then printf '%s\n' "$3";else printf '%s\n' "$v";fi; }
"$4" "$5"
'''
        env = dict(os.environ, PATH=str(Path(sys.executable).parent) + os.pathsep + os.environ["PATH"])
        return subprocess.run(["bash", "--noprofile", "--norc", "-c", body, "fixture",
            str(ROOT / "scripts/install_utils.sh"), str(self.scripts), str(self.folder), function, str(self.config)],
            env=env, text=True, capture_output=True, timeout=30)

    def test_both_real_wrappers_preserve_options_and_completed_warm_reuse(self):
        for kind in M.PRODUCERS:
            with self.subTest(kind=kind):
                args = self.setup_wrapper(kind)
                result = self.wrapper(kind)
                self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
                call = json.loads(self.argv_log.read_text().splitlines()[-1])
                self.assertEqual(call[call.index("--nm_to_enst_tsv") + 1], str(self.mapping))
                self.assertEqual("--skip_das" in call, kind == "pm1")
                self.assertEqual(call[call.index("--threads" if kind == "pm1" else "--processes") + 1], "2")
                self.assertEqual(pickle.loads(args.artifact.read_bytes()), PAYLOAD)
                before = self.pair(args.artifact), self.argv_log.read_bytes()
                result = self.wrapper(kind)
                self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
                self.assertEqual((self.pair(args.artifact), self.argv_log.read_bytes()), before)

    def test_both_wrappers_recover_before_receipt_gate_and_do_not_recalculate(self):
        for kind in M.PRODUCERS:
            with self.subTest(kind=kind):
                args = self.setup_wrapper(kind)
                self.interrupt(args, "artifact")
                before = self.argv_log.read_bytes()
                result = self.wrapper(kind)
                self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
                self.assertEqual(self.argv_log.read_bytes(), before)
                M.BIOMART.check_receipt(args.artifact, self.mapping, self.anchors)

    def test_both_wrappers_refuse_unknown_existing_artifact_without_config_change(self):
        for kind in M.PRODUCERS:
            with self.subTest(kind=kind):
                args = self.setup_wrapper(kind)
                args.artifact.write_bytes(b"unknown supplied artifact")
                before = self.config.read_bytes(), M.file_state(args.artifact)
                result = self.wrapper(kind)
                self.assertNotEqual(result.returncode, 0)
                self.assertIn("missing/stale NM-to-ENST receipt", result.stderr)
                self.assertEqual((self.config.read_bytes(), M.file_state(args.artifact)), before)
                self.assertFalse(self.argv_log.exists())

    def test_pair_committed_before_failed_yaml_update_reuses_without_calculation(self):
        for kind in M.PRODUCERS:
            with self.subTest(kind=kind):
                args = self.setup_wrapper(kind)
                before_config = self.config.read_bytes()
                with mock.patch.dict(os.environ, {"FAIL_CONFIG_WRITE": "1"}):
                    result = self.wrapper(kind)
                self.assertNotEqual(result.returncode, 0)
                self.assertEqual(self.config.read_bytes(), before_config)
                M.BIOMART.check_receipt(args.artifact, self.mapping, self.anchors)
                before = self.pair(args.artifact), self.argv_log.read_bytes()
                result = self.wrapper(kind)
                self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
                self.assertEqual((self.pair(args.artifact), self.argv_log.read_bytes()), before)


if __name__ == "__main__":
    unittest.main()
