"""Private builder-view routing with no annotation, downloads or full installation.

The shell boundary uses a labelled bootstrap and cache-install double. An
optional existing compiled CLI exercises the real prepare/publish bridge too.
"""
import hashlib
import json
import os
from pathlib import Path
import pickle
import subprocess
import sys
import tempfile
import time
import unittest

import yaml


ROOT = Path(__file__).resolve().parents[1]
INSTALLER = ROOT / "scripts/install_utils.sh"
WRAPPED = (
    "AlphaMissense_stat", "AlphaMissense_pick_missense_constrained_segments",
    "AlphaMissense_pick_intolerant_domains", "prepare_pm1_regions", "InterPro_parsing",
    "ClinGen_deploy", "ClinVar_AA_stat", "ClinVar_patho_AF_stat", "ClinVar_Gene_stat",
    "ClinVar_Gene_stat_regenerate", "ClinVar_install",
)
STAGES = (
    "priva_manual_cache_inputs", "priva_prepare_install_config", "bundled_snapshot_install",
    "Nextflow_runtime_install", "PriVA_rust_install", "reference_genome_install",
    "priva_verify_vep_release", "vep_cache_api_install", "VEP_plugins_install", "LoFtee_install",
    "gnomAD_install", "AlphaMissense_anno_gnomAD", "InterPro_parsing", "pext_install",
    "ClinVar_install", "AlphaMissense_stat", "AlphaMissense_pick_intolerant_domains",
    "prepare_pm1_regions", "ClinVar_patho_AF_stat", "ClinVar_PM5_index", "CADD_install",
    "CADD_upstream_check", "ClinGen_deploy", "CDS_FASTA_install", "mechanism_resource_install",
    "priva_export_installed_environment",
)
FAKE_BINARY = r'''#!/usr/bin/env python3
import hashlib,json,os,sys
from pathlib import Path
import yaml
args=sys.argv[1:]
def event(*items):
    with open(os.environ['EVENTS'],'a') as output: output.write(json.dumps(items)+'\n')
if args == ['--version']:
    print('priva fixture'); raise SystemExit(0)
assert args[0] == 'cache-install'
operation=args[1]
event(operation,*args[2:])
if os.environ.get('FAIL_CACHE') == operation: raise SystemExit(61)
values=dict(zip(args[2::2],args[3::2]))
original=Path(values['--config']); private=Path(values['--build-config'])
context=Path(str(private)+'.cache-install.json')
sha=lambda: hashlib.sha256(original.read_bytes()).hexdigest()
if operation == 'prepare':
    if context.exists():
        prior=json.loads(context.read_text())
        if prior['state']=='pending':
            if prior['sha256'] != sha(): raise SystemExit(62)
            raise SystemExit(0)
    config=yaml.safe_load(original.read_text())
    for key in ('all_intolerant_domains','clinvar_aa_stat'):
        if isinstance(config.get(key),str) and config[key].endswith('.rkyv'):
            config[key]=config[key][:-5]+'.pkl'
    private.write_text(yaml.safe_dump(config,sort_keys=False))
    context.write_text(json.dumps({'sha256':sha(),'state':'pending'}))
elif operation == 'publish':
    prior=json.loads(context.read_text())
    if prior['sha256'] != sha(): raise SystemExit(63)
    config=yaml.safe_load(private.read_text())
    for key in ('all_intolerant_domains','clinvar_aa_stat'):
        if isinstance(config.get(key),str) and config[key].endswith('.pkl'):
            config[key]=config[key][:-4]+'.rkyv'
    temporary=original.with_name(original.name+'.native.tmp')
    temporary.write_text(yaml.safe_dump(config,sort_keys=False)); temporary.replace(original)
    context.write_text(json.dumps({'sha256':sha(),'state':'complete'}))
else: raise AssertionError(operation)
'''
PRELUDE = r'''
source "$INSTALLER"
wget() { echo UNEXPECTED_DOWNLOAD >&2; return 98; }
curl() { echo UNEXPECTED_DOWNLOAD >&2; return 98; }
record_event() {
    "$TEST_PYTHON" -c 'import json,os,sys;open(os.environ["EVENTS"],"a").write(json.dumps(sys.argv[1:])+"\n")' "$@"
}
read_yaml() {
    "$TEST_PYTHON" -c 'import sys,yaml;v=yaml.safe_load(open(sys.argv[1])).get(sys.argv[2]);print("null" if v is None else v)' "$1" "$2"
}
python() {
    if [[ $1 == */bootstrap_priva_rust.py ]]; then
        record_event bootstrap "$@"
        [[ ${FAIL_CACHE:-} != bootstrap ]] || return 60
        "$TEST_PYTHON" -c 'import json,os;print(json.dumps({"binary":os.environ["FAKE_BIN"]}))'
    else "$TEST_PYTHON" "$@"; fi
}
fixture_builder() {
    record_event builder "$@"
    [[ ${FAIL_BUILDER:-0} != 1 ]] || return 64
}
'''


class NativeCacheInstallWiring(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory(prefix="native_install_wiring_",
            dir=os.environ.get("PRIVA_TEST_TMPDIR", os.environ.get("TMPDIR")))
        self.addCleanup(self.temporary.cleanup)
        self.root = Path(self.temporary.name)
        self.cache = self.root / "cache"
        self.cache.mkdir()
        self.tmp = self.root / "tmp"
        self.tmp.mkdir()
        self.events = self.root / "events.jsonl"
        self.events.write_text("")
        self.binary = self.root / "priva-fixture"
        self.binary.write_text(FAKE_BINARY.replace("#!/usr/bin/env python3", "#!" + sys.executable, 1))
        self.binary.chmod(0o755)
        self.config = self.root / "runtime config.yaml"
        self.values = {"cache_root": str(self.cache), "tmp_dir": str(self.tmp),
            "assembly": "hg38", "vep_installed_version": 116,
            "all_intolerant_domains": str(self.cache / "domains.v1.rkyv"),
            "priva_binary": "previous-selected-program"}
        self.write_config()
        self.env = dict(os.environ, INSTALLER=str(INSTALLER), TEST_PYTHON=sys.executable,
                        EVENTS=str(self.events), FAKE_BIN=str(self.binary), CONFIG=str(self.config))

    def write_config(self):
        self.config.write_text(yaml.safe_dump(self.values, sort_keys=False))

    def execute(self, body, **extra):
        result = subprocess.run(["bash", "--noprofile", "--norc", "-c", PRELUDE + body],
            env=dict(self.env, **extra), cwd=self.root, capture_output=True, text=True, timeout=40)
        self.assertNotIn("UNEXPECTED_DOWNLOAD", result.stdout + result.stderr)
        events = [json.loads(line) for line in self.events.read_text().splitlines()]
        return result, events

    def private_path(self):
        identity = hashlib.sha256(str(self.config.resolve()).encode()).hexdigest()
        return self.cache / "runtime/priva-cache-configs" / identity / "build.yaml"

    def test_legacy_and_unrelated_native_names_keep_original_arguments(self):
        self.values["all_intolerant_domains"] = str(self.cache / "selected.pkl")
        self.values["notes"] = "A comment about .rkyv is not a configured native cache"
        self.write_config()
        result, events = self.execute(r'''
_priva_cache_builder_ClinVar_Gene_stat_regenerate() { fixture_builder "$@"; }
ClinVar_Gene_stat_regenerate "$CONFIG" "summary with spaces.gz" mane.gz "candidate's.pkl.gz"
''')
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertEqual(events, [["builder", str(self.config), "summary with spaces.gz", "mane.gz", "candidate's.pkl.gz"]])

    def test_all_direct_builder_entrypoints_receive_private_pickle_view(self):
        for function in WRAPPED:
            with self.subTest(function=function):
                self.events.write_text("")
                self.write_config()
                position = 1 if function == "ClinVar_install" else 0
                call = f'{function} --config "$CONFIG"' if position else f'{function} "$CONFIG" "extra argument"'
                body = f'_priva_cache_builder_{function}() {{ fixture_builder "$@"; }}\n{call}\n'
                result, events = self.execute(body)
                self.assertEqual(result.returncode, 0, result.stderr)
                self.assertEqual([row[0] for row in events], ["bootstrap", "prepare", "builder", "publish"])
                args = events[2][1:]
                self.assertEqual(args[position], str(self.private_path()))
                self.assertTrue(yaml.safe_load(self.private_path().read_text())["all_intolerant_domains"].endswith(".pkl"))
                self.assertTrue(yaml.safe_load(self.config.read_text())["all_intolerant_domains"].endswith(".rkyv"))

    def test_nested_private_builder_does_not_prepare_or_publish_again(self):
        result, events = self.execute(r'''
_priva_cache_builder_ClinVar_AA_stat() { fixture_builder "$@"; }
_priva_cache_builder_ClinVar_install() { ClinVar_AA_stat "$2"; }
ClinVar_install --config "$CONFIG"
''')
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertEqual([row[0] for row in events], ["bootstrap", "prepare", "builder", "publish"])

    def test_repeated_wrapper_registration_retains_the_raw_builder(self):
        result, events = self.execute(r'''
_priva_cache_builder_ClinVar_AA_stat() { fixture_builder "$@"; }
priva_wrap_cache_builders
ClinVar_AA_stat "$CONFIG"
''')
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertEqual([row[0] for row in events], ["bootstrap", "prepare", "builder", "publish"])

    def test_prepare_builder_and_publish_failures_preserve_original_config(self):
        for failure in ("bootstrap", "prepare", "builder", "publish"):
            with self.subTest(failure=failure):
                self.events.write_text("")
                self.write_config()
                before = self.config.read_bytes()
                result, events = self.execute(r'''
_priva_cache_builder_ClinVar_AA_stat() { fixture_builder "$@"; }
ClinVar_AA_stat "$CONFIG"
''', FAIL_CACHE=failure, FAIL_BUILDER="1" if failure == "builder" else "0")
                self.assertNotEqual(result.returncode, 0)
                self.assertEqual(self.config.read_bytes(), before)
                if failure in ("bootstrap", "prepare"):
                    self.assertNotIn("builder", [row[0] for row in events])
                if failure == "builder":
                    self.assertNotIn("publish", [row[0] for row in events])

    def test_failed_builder_resumes_same_private_file_and_releases_lock(self):
        body = r'''
_priva_cache_builder_ClinVar_AA_stat() {
    "$TEST_PYTHON" - "$1" <<'PY'
from pathlib import Path
import sys,yaml
p=Path(sys.argv[1]);v=yaml.safe_load(p.read_text());v['attempt']=v.get('attempt',0)+1;p.write_text(yaml.safe_dump(v))
PY
    fixture_builder "$@"
}
ClinVar_AA_stat "$CONFIG"
'''
        before = self.config.read_bytes()
        result, _ = self.execute(body, FAIL_BUILDER="1")
        self.assertNotEqual(result.returncode, 0)
        self.assertEqual(self.config.read_bytes(), before)
        self.assertEqual(yaml.safe_load(self.private_path().read_text())["attempt"], 1)
        result, _ = self.execute(body, FAIL_BUILDER="0")
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertEqual(yaml.safe_load(self.config.read_text())["attempt"], 2)

    def test_second_caller_cannot_prepare_while_first_builder_holds_lock(self):
        ready = self.root / "first-builder-ready"
        release = self.root / "release-first-builder"
        body = r'''
_priva_cache_builder_ClinVar_AA_stat() {
    printf '\nprivate_progress: retained\n' >> "$1"
    record_event builder_entered "$@"
    : > "$READY"
    while [[ ! -e $RELEASE ]]; do sleep 0.05; done
    fixture_builder "$@"
}
ClinVar_AA_stat "$CONFIG"
'''
        first = subprocess.Popen(["bash", "--noprofile", "--norc", "-c", PRELUDE + body],
            env=dict(self.env, READY=str(ready), RELEASE=str(release)), cwd=self.root,
            stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
        original = self.config.read_bytes()
        try:
            deadline = time.monotonic() + 15
            while not ready.exists() and first.poll() is None and time.monotonic() < deadline:
                time.sleep(0.05)
            self.assertTrue(ready.exists(), "first builder did not reach its held-lock phase")
            private = self.private_path().read_bytes()
            context = Path(str(self.private_path()) + ".cache-install.json").read_bytes()
            second, events = self.execute(r'''
_priva_cache_builder_ClinVar_AA_stat() { fixture_builder "$@"; }
ClinVar_AA_stat "$CONFIG"
''')
            self.assertNotEqual(second.returncode, 0)
            self.assertIn("Another installation", second.stderr)
            self.assertEqual(sum(row[0] == "prepare" for row in events), 1)
            self.assertEqual(self.private_path().read_bytes(), private)
            self.assertEqual(Path(str(self.private_path()) + ".cache-install.json").read_bytes(), context)
            self.assertEqual(self.config.read_bytes(), original)
        finally:
            release.write_text("release\n")
            try:
                stdout, stderr = first.communicate(timeout=15)
            except subprocess.TimeoutExpired:
                first.kill()
                first.communicate()
                raise
        self.assertEqual(first.returncode, 0, stdout + stderr)

    def test_user_edit_during_builder_is_not_overwritten(self):
        result, events = self.execute(r'''
_priva_cache_builder_ClinVar_AA_stat() {
    printf '\nuser_edit: retained\n' >> "$CONFIG"
    fixture_builder "$@"
}
ClinVar_AA_stat "$CONFIG"
''')
        self.assertNotEqual(result.returncode, 0)
        self.assertEqual(yaml.safe_load(self.config.read_text())["user_edit"], "retained")
        self.assertEqual([row[0] for row in events], ["bootstrap", "prepare", "builder", "publish"])

    def test_private_parent_symlink_is_rejected(self):
        outside = self.root / "unselected"
        outside.mkdir()
        (self.cache / "runtime").symlink_to(outside, target_is_directory=True)
        result, events = self.execute(r'''
_priva_cache_builder_ClinVar_AA_stat() { fixture_builder "$@"; }
ClinVar_AA_stat "$CONFIG"
''')
        self.assertNotEqual(result.returncode, 0)
        self.assertIn("real managed directory", result.stderr)
        self.assertEqual(list(outside.iterdir()), [])
        self.assertNotIn("prepare", [row[0] for row in events])

    def test_private_lock_cannot_truncate_a_linked_file(self):
        self.private_path().parent.mkdir(parents=True)
        lock = Path(str(self.private_path()) + ".lock")
        lock.symlink_to(self.config)
        before = self.config.read_bytes()
        result, events = self.execute(r'''
_priva_cache_builder_ClinVar_AA_stat() { fixture_builder "$@"; }
ClinVar_AA_stat "$CONFIG"
''')
        self.assertNotEqual(result.returncode, 0)
        self.assertEqual(self.config.read_bytes(), before)
        self.assertNotIn("prepare", [row[0] for row in events])

    def test_main_uses_private_config_and_publishes_only_after_all_stages(self):
        for key in ("vep_cache_dir", "vep_plugins_dir", "vep_plugins_cachedir"):
            path = self.cache / key
            path.mkdir()
            self.values[key] = str(path)
        reference = self.root / "reference.fa"
        reference.write_text(">chr1\nA\n")
        self.values["ref_genome"] = str(reference)
        self.write_config()
        body = 'conda_install_vep() { PRIVA_INSTALL_ENV_NAME=priva_acmg; }\n'
        for stage in STAGES:
            target = "_priva_cache_builder_" + stage if stage in WRAPPED else stage
            body += f'{target}() {{ record_event {stage} "$@"; }}\n'
        body += 'main_install "$CONFIG"\n'
        result, events = self.execute(body)
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertEqual([row[0] for row in events], ["bootstrap", "prepare", *STAGES, "publish"])
        for row in events[2:-1]:
            if row[0] in ("priva_verify_vep_release", "vep_cache_api_install"):
                continue
            self.assertIn(str(self.private_path()), row[1:], row[0])

    def test_main_builder_failure_retains_private_progress_without_publication(self):
        for key in ("vep_cache_dir", "vep_plugins_dir", "vep_plugins_cachedir"):
            path = self.cache / key
            path.mkdir()
            self.values[key] = str(path)
        reference = self.root / "reference.fa"
        reference.write_text(">chr1\nA\n")
        self.values["ref_genome"] = str(reference)
        self.write_config()
        before = self.config.read_bytes()
        body = 'conda_install_vep() { PRIVA_INSTALL_ENV_NAME=priva_acmg; }\n'
        for stage in STAGES:
            target = "_priva_cache_builder_" + stage if stage in WRAPPED else stage
            body += f'{target}() {{ record_event {stage} "$@"; }}\n'
        body += r'''
_priva_cache_builder_AlphaMissense_stat() {
    printf '\ncompleted_private_stage: retained\n' >> "$1"
    record_event failed_builder "$@"
    return 65
}
main_install "$CONFIG"
'''
        result, events = self.execute(body)
        self.assertNotEqual(result.returncode, 0)
        self.assertEqual(self.config.read_bytes(), before)
        self.assertEqual(yaml.safe_load(self.private_path().read_text())["completed_private_stage"], "retained")
        self.assertNotIn("publish", [row[0] for row in events])

    @unittest.skipUnless(os.environ.get("PRIVA_CACHE_INSTALL_TEST_BINARY"), "an existing native CLI is optional")
    def test_real_prepare_publish_bridge_preserves_native_archive_on_warm_builder(self):
        binary = Path(os.environ["PRIVA_CACHE_INSTALL_TEST_BINARY"])
        source = self.cache / "source.pkl"
        source.write_bytes(pickle.dumps({"domain:a", "domain:b"}, protocol=4))
        self.values["all_intolerant_domains"] = str(source)
        self.write_config()
        private = self.root / "initial builder.yaml"
        for operation in ("prepare", "publish"):
            subprocess.run([str(binary), "cache-install", operation, "--config", str(self.config),
                "--build-config", str(private), "--source-root", str(ROOT)], check=True, capture_output=True, text=True)
        archive = Path(yaml.safe_load(self.config.read_text())["all_intolerant_domains"])
        before = archive.read_bytes(), archive.stat().st_mtime_ns
        self.env["FAKE_BIN"] = str(binary)
        result, events = self.execute(r'''
_priva_cache_builder_ClinVar_AA_stat() { fixture_builder "$@"; }
ClinVar_AA_stat "$CONFIG"
''')
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertEqual((archive.read_bytes(), archive.stat().st_mtime_ns), before)
        self.assertEqual([row[0] for row in events], ["bootstrap", "builder"])
        self.assertEqual(yaml.safe_load(self.private_path().read_text())["all_intolerant_domains"], str(source))


if __name__ == "__main__":
    unittest.main()
