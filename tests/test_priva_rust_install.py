"""Managed CLI publication/reuse controls; no real download, Cargo or models.

Tests exercise real byte copying, qualified-receipt binding, no-replace rename,
read-only reuse, source/evidence races and the actual shell install function.
The separate deployment check imports the already qualified real executable.
"""
import importlib.util
import io
import json
import os
from pathlib import Path
import subprocess
import sys
import tempfile
from types import SimpleNamespace
import unittest
from unittest import mock

ROOT = Path(__file__).resolve().parents[1]
SPEC = importlib.util.spec_from_file_location("managed_priva_rust", ROOT / "scripts/bootstrap_priva_rust.py")
B = importlib.util.module_from_spec(SPEC); SPEC.loader.exec_module(B)


def write_json(path, value):
    path.write_bytes(B.canonical(value) + b"\n")


def elf(extra=b"fixture"):
    header = bytearray(64); header[:6] = b"\x7fELF\x02\x01"; header[18:20] = b"\x3e\x00"
    return bytes(header) + extra


class ManagedRustTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory(dir=os.environ.get("PRIVA_TEST_TMPDIR", os.environ.get("TMPDIR")))
        self.addCleanup(self.temp.cleanup); self.folder = Path(self.temp.name)
        self.source = self.folder / "repo"; (self.source / "rust/src").mkdir(parents=True)
        (self.source / "rust/src/main.rs").write_text("fn main() {}\n")
        (self.source / "rust/src/PANDAS_LICENSE").write_text("owned source-license fixture\n")
        (self.source / "rust/src/EXPAT_LICENSE").write_text("static XML dependency notice fixture\n")
        (self.source / "rust/src/dependency").mkdir()
        (self.source / "rust/src/dependency/PANDAS_LICENSE").write_text("distinct nested notice fixture\n")
        (self.source / "rust/third_party").mkdir()
        (self.source / "rust/third_party/parquet_licenses.txt").write_text("Parquet dependency notice fixture\n")
        (self.source / "rust/Cargo.toml").write_text('''[package]
name = "priva-tabular"
version = "0.1.0"
[[bin]]
name = "priva"
path = "src/main.rs"
[profile.dev]
debug = 0
[profile.test]
debug = 0
''')
        (self.source / "rust/Cargo.lock").write_text("version = 4\n")
        (self.source / "rust/rust-toolchain.toml").write_text('[toolchain]\nchannel = "1.90.0"\nprofile = "minimal"\n')
        (self.source / "data/liftover").mkdir(parents=True)
        (self.source / "data/liftover/GRC_to_ucsc.contig.map.txt").write_text("1 chr1\n")
        (self.source / "scripts").mkdir()
        (self.source / "scripts/das_input_continuity.py").write_text("# embedded checker fixture\n")
        for name in ("biomart_nm_to_enst.py", "biomart_ftp_nm_to_enst.py"):
            (self.source / "scripts" / name).write_text("# embedded BioMart source fixture\n")
        self.cache = self.folder / "cache"; self.cache.mkdir()
        self.binary = self.folder / "qualified-priva"; self.binary.write_bytes(elf()); self.binary.chmod(0o700)
        self.pins_path = self.folder / "pins.json"; self.validation_path = self.folder / "validation.json"
        self.pins = dict(source_files=B.source_state(self.source)["files"], actual_Cargo_profile=B.IMPORTED_PROFILE,
            toolchain={"/frozen/1.90.0/bin/" + name: B.sha(name.encode()) for name in ("cargo", "rustc", "rustdoc", "rustfmt")})
        self.result = dict(completed=True, source_pins_sha256="", rust_tests_passed=162,
            binary_sha256=B.digest(self.binary), native=dict(job_id="closed-fixture", step_id="0"),
            scope="explicit fixture procedure qualification", scientific_acceptance=False, production_backend_selected=False)
        self.refresh()
        self.probe = mock.patch.object(B, "probe", return_value="priva 0.1.0 (fixture)").start()
        self.addCleanup(mock.patch.stopall)

    def refresh(self):
        write_json(self.pins_path, self.pins)
        self.result["source_pins_sha256"] = B.digest(self.pins_path)
        write_json(self.validation_path, self.result)

    def args(self, imported=True, **changes):
        result = dict(source_root=self.source, cache_root=self.cache, identity=None,
            qualified_binary=self.binary if imported else None, build_pins=self.pins_path if imported else None,
            build_pins_sha256=B.digest(self.pins_path) if imported else None,
            validation=self.validation_path if imported else None,
            validation_sha256=B.digest(self.validation_path) if imported else None,
            toolchain_bin=None, build_profile="release", jobs=1, offline=False)
        result.update(changes); return SimpleNamespace(**result)

    def snapshot(self, root):
        return {str(path.relative_to(root)): (B.identity(path), path.read_bytes())
                for path in root.rglob("*") if path.is_file() and not path.is_symlink()}

    def test_real_copy_import_retains_profile_scope_and_source_identity(self):
        before = B.identity(self.binary); result = B.install(self.args())
        binary = Path(result["binary"])
        self.assertEqual(binary.read_bytes(), self.binary.read_bytes())
        self.assertNotEqual(binary.stat().st_ino, self.binary.stat().st_ino)
        self.assertEqual(B.identity(self.binary), before)
        self.assertEqual(result["profile"]["name"], "dev")
        self.assertEqual(result["profile"]["opt_level"], 0)
        self.assertFalse(result["qualification"]["scientific_acceptance"])
        self.assertFalse(result["qualification"]["production_backend_selected"])
        receipt, _, _ = B.load(Path(result["receipt"]))
        self.assertEqual(set(receipt["notices_sha256"]), {
            "licenses/PANDAS_LICENSE", "licenses/EXPAT_LICENSE", "licenses/dependency/PANDAS_LICENSE",
            "licenses/third_party/parquet_licenses.txt"})
        for relative, wanted in receipt["notices_sha256"].items():
            installed = binary.parent.parent / relative
            original = self.source / B.source_notice_paths(B.source_state(self.source)['files'])[relative]
            self.assertEqual(installed.read_bytes(), original.read_bytes())
            self.assertEqual(B.digest(installed), wanted)
            self.assertNotEqual(installed.stat().st_ino, original.stat().st_ino)

    def test_warm_install_and_check_are_read_only_and_do_not_build(self):
        result = B.install(self.args()); before = self.snapshot(self.cache)
        with mock.patch.object(B, "cold_build", side_effect=AssertionError("unexpected build")), \
             mock.patch.object(B, "obtain_toolchain", side_effect=AssertionError("unexpected toolchain")):
            self.assertEqual(B.install(self.args()), result)
            self.assertEqual(B.install(self.args(imported=False)), result)
            self.assertEqual(B.check(self.args(imported=False)), result)
        self.assertEqual(self.snapshot(self.cache), before)

    def test_check_absent_does_not_create_runtime(self):
        with self.assertRaises(ValueError): B.check(self.args(imported=False))
        self.assertEqual(list(self.cache.iterdir()), [])

    def test_explicit_receipt_sha_change_refused(self):
        with self.assertRaises(ValueError): B.install(self.args(build_pins_sha256="0" * 64))
        self.assertEqual(list(self.cache.iterdir()), [])

    def test_partial_import_refused(self):
        with self.assertRaises(ValueError): B.install(self.args(validation=None))

    def test_unfinished_validation_refused(self):
        self.result["completed"] = False; self.refresh()
        with self.assertRaises(ValueError): B.install(self.args())

    def test_validation_must_bind_exact_build_pins(self):
        self.result["source_pins_sha256"] = "0" * 64; write_json(self.validation_path, self.result)
        with self.assertRaises(ValueError): B.install(self.args())

    def test_different_cargo_lock_and_extra_rust_source_refused(self):
        for name, raw in (("rust/Cargo.lock", "version = 3\n"), ("rust/src/extra.rs", "pub fn extra() {}\n")):
            with self.subTest(name=name):
                path = self.source / name; before = path.read_bytes() if path.exists() else None
                path.write_text(raw)
                with self.assertRaises(ValueError): B.install(self.args())
                if before is None: path.unlink()
                else: path.write_bytes(before)

    def test_non_rs_source_license_is_pinned_and_changes_invalidate_import(self):
        state = B.source_state(self.source)
        self.assertIn('rust/src/PANDAS_LICENSE', state['files'])
        self.assertIn('rust/third_party/parquet_licenses.txt', state['files'])
        for name in ('rust/src/PANDAS_LICENSE', 'rust/third_party/parquet_licenses.txt'):
            with self.subTest(name=name):
                path = self.source / name
                original = path.read_bytes()
                path.write_text('changed source companion\n')
                with self.assertRaises(ValueError): B.install(self.args())
                path.write_bytes(original)

    def test_new_build_profile_cannot_be_relabelled_as_qualified(self):
        self.pins["actual_Cargo_profile"] = "release"; self.refresh()
        with self.assertRaises(ValueError): B.install(self.args())

    def test_changed_binary_blocks_publication_and_retains_failure(self):
        self.binary.write_bytes(elf(b"changed"))
        with self.assertRaises(ValueError): B.install(self.args())
        base = self.cache / "runtime/priva-rust"
        self.assertFalse(any(len(path.name) == 64 for path in base.iterdir()))
        self.assertEqual(len(list(base.glob(".artifact-work-*/failure.json"))), 1)

    def test_non_elf_and_literal_import_links_refused(self):
        self.binary.write_bytes(b"#!/bin/sh\nexit 0\n"); self.result["binary_sha256"] = B.digest(self.binary); self.refresh()
        with self.assertRaises(ValueError): B.install(self.args())
        link = self.folder / "binary-link"; link.symlink_to(self.binary.name)
        with self.assertRaises(ValueError): B.install(self.args(qualified_binary=link))

    def test_proof_race_during_actual_copy_and_probe_refused(self):
        def change(*_):
            self.validation_path.write_bytes(self.validation_path.read_bytes() + b" ")
            return "priva 0.1.0 (fixture)"
        self.probe.side_effect = change
        with self.assertRaises(ValueError): B.install(self.args())
        self.assertFalse(any(len(path.name) == 64 for path in (self.cache / "runtime/priva-rust").iterdir()))

    def test_source_race_before_publication_refused(self):
        def change(*_):
            (self.source / "rust/src/main.rs").write_text("fn main() { panic!(); }\n")
            return "priva 0.1.0 (fixture)"
        self.probe.side_effect = change
        with self.assertRaises(ValueError): B.install(self.args())
        self.assertFalse(any(len(path.name) == 64 for path in (self.cache / "runtime/priva-rust").iterdir()))

    def test_probe_failure_never_marks_partial_ready(self):
        self.probe.side_effect = subprocess.CalledProcessError(2, ["fixture", "--version"])
        with self.assertRaises(subprocess.CalledProcessError): B.install(self.args())
        self.assertEqual(len(list((self.cache / "runtime/priva-rust").glob(".artifact-work-*/failure.json"))), 1)

    def test_final_publication_never_clobbers_a_racing_target(self):
        publish = B.publish_artifact
        def race(source, target):
            target.mkdir(); (target / "foreign").write_text("preserve")
            return publish(source, target)
        with mock.patch.object(B, "publish_artifact", side_effect=race):
            with self.assertRaises(FileExistsError): B.install(self.args())
        self.assertEqual([p.read_text() for p in (self.cache / "runtime/priva-rust").glob("*/foreign")], ["preserve"])

    def test_empty_racing_target_never_adopted_by_either_retry_route(self):
        for explicit in (True, False):
            with self.subTest(explicit_import=explicit):
                self.cache = self.folder / ('empty-race-' + str(explicit)); self.cache.mkdir()
                publish = B.publish_artifact; foreign = []
                def race(stage, target):
                    target.mkdir(mode=0o700); foreign.append(target)
                    return publish(stage, target)
                with mock.patch.object(B, 'publish_artifact', side_effect=race):
                    with self.assertRaises(FileExistsError): B.install(self.args())
                target = foreign[0]; before = B.identity(target)
                pending_before = self.snapshot(self.cache)
                self.assertEqual(list(target.iterdir()), [])
                result = None; rejected = None
                try: result = B.install(self.args(imported=explicit))
                except ValueError as error: rejected = error
                self.assertIsNone(result, 'retry adopted an empty foreign target: ' + repr(result))
                self.assertIsNotNone(rejected)
                self.assertEqual(B.identity(target), before)
                self.assertEqual(list(target.iterdir()), [])
                self.assertEqual(self.snapshot(self.cache), pending_before)

    def test_completed_old_publication_remains_read_only_reusable(self):
        ready = B.install(self.args()); target = Path(ready['receipt']).parent
        ledger_path = B.publication_path(target); ledger, _, _ = B.load(ledger_path)
        # The first deployed protocol has a complete byte-bound ready receipt
        # but predates notices and target-creation observations. Fixture only.
        receipt, _, _ = B.load(target / 'receipt.json')
        notices = receipt.pop('notices_sha256')
        for name in notices:
            (target / name).unlink(); ledger['members'].pop(name)
        for folder in sorted(B.member_directories(notices), key=lambda name: len(Path(name).parts), reverse=True):
            (target / folder).rmdir()
        write_json(target / 'receipt.json', receipt)
        ledger['members']['receipt.json']['sha256'] = B.digest(target / 'receipt.json')
        ledger.pop('protocol_version', None); write_json(ledger_path, ledger)
        created = target.parent / ('publication-' + target.name + '.created.json')
        if created.exists(): created.unlink()
        before = self.snapshot(self.cache)
        self.assertEqual(B.install(self.args()), ready)
        self.assertEqual(B.check(self.args(imported=False)), ready)
        self.assertEqual(self.snapshot(self.cache), before)

    def test_interruption_before_creation_proof_fails_closed_on_retry(self):
        atomic_new = B.atomic_new
        def interrupt(path, raw):
            if str(path).endswith('.created.json'):
                raise OSError('fixture interruption before target ownership was committed')
            return atomic_new(path, raw)
        with mock.patch.object(B, 'atomic_new', side_effect=interrupt):
            with self.assertRaises(OSError): B.install(self.args())
        target = next(path for path in (self.cache/'runtime/priva-rust').iterdir() if len(path.name) == 64)
        before = B.identity(target), self.snapshot(self.cache)
        with self.assertRaises(ValueError): B.install(self.args())
        self.assertEqual((B.identity(target), self.snapshot(self.cache)), before)
        self.assertEqual(list(target.iterdir()), [])

    def test_replaced_owned_directory_cannot_resume_into_empty_foreign_target(self):
        target = self.interrupt_publication()
        retained = target.parent / '.retained-own-directory'; target.rename(retained)
        target.mkdir(mode=0o700)
        before = B.identity(target), self.snapshot(self.cache)
        with self.assertRaises(ValueError): B.install(self.args())
        self.assertEqual((B.identity(target), self.snapshot(self.cache)), before)
        self.assertEqual(list(target.iterdir()), [])

    def test_creation_proof_ledger_binding_is_required(self):
        target = self.interrupt_publication()
        proof = target.parent / ('publication-' + target.name + '.created.json')
        value, _, _ = B.load(proof); value['publication_sha256'] = '0' * 64
        write_json(proof, value); before = self.snapshot(self.cache)
        with self.assertRaises(ValueError): B.install(self.args())
        self.assertEqual(self.snapshot(self.cache), before)
        self.assertFalse((target/'receipt.json').exists())

    def interrupt_publication(self):
        link = os.link
        def interrupted(source, target, **kwargs):
            if Path(target).parent.name == 'evidence':
                raise OSError('fixture interruption after the binary was linked')
            return link(source, target, **kwargs)
        with mock.patch.object(B.os, 'link', side_effect=interrupted):
            with self.assertRaises(OSError): B.install(self.args())
        base = self.cache / 'runtime/priva-rust'
        target = next(path for path in base.iterdir() if len(path.name) == 64)
        self.assertFalse((target/'receipt.json').exists())
        return target

    def test_interrupted_owned_publication_recovers_exactly_without_recopy(self):
        target = self.interrupt_publication()
        before = self.snapshot(self.cache)
        with self.assertRaises(ValueError): B.check(self.args(imported=False))
        self.assertEqual(self.snapshot(self.cache), before)
        with mock.patch.object(B, 'copy_binary', side_effect=AssertionError('must reuse prepared bytes')):
            ready = B.install(self.args())
        self.assertEqual(ready['identity'], target.name)
        self.assertTrue((target/'receipt.json').is_file())
        self.assertEqual(len(list((self.cache/'runtime/priva-rust').glob('.artifact-work-*/failure.json'))), 1)

    def test_changed_owned_partial_is_refused_and_source_elf_preserved(self):
        target = self.interrupt_publication(); before = B.identity(self.binary)
        path = target/'bin/priva'; path.chmod(0o755); path.write_bytes(elf(b'changed partial'))
        with self.assertRaises(ValueError): B.install(self.args())
        self.assertFalse((target/'receipt.json').exists())
        self.assertEqual(B.identity(self.binary), before)

    def test_existing_root_alias_supported_and_child_alias_refused(self):
        alias = self.folder / "cache-alias"; alias.symlink_to(self.cache.name, target_is_directory=True)
        result = B.install(self.args(cache_root=alias)); self.assertTrue(result["ready"])
        other = self.folder / "other"; other.mkdir()
        cache = self.folder / "second-cache"; cache.mkdir(); (cache / "runtime").symlink_to(other)
        with self.assertRaises(ValueError): B.install(self.args(cache_root=cache))
        self.assertEqual(list(other.iterdir()), [])

    def test_foreign_runtime_and_malformed_receipt_are_retained(self):
        base = self.cache / "runtime/priva-rust"; base.mkdir(parents=True)
        (base / "foreign").write_bytes(b"retained")
        with self.assertRaises((ValueError, OSError)): B.install(self.args())
        self.assertEqual((base / "foreign").read_bytes(), b"retained")

    def test_modified_executable_evidence_and_extra_members_fail_warm_check(self):
        ready = B.install(self.args()); root = Path(ready["receipt"]).parent
        for relative, raw in (("bin/priva", elf(b"different")), ("evidence/validation.json", b"{}\n"),
                              ("licenses/EXPAT_LICENSE", b"changed notice"), ("foreign", b"unexpected")):
            with self.subTest(relative=relative):
                path = root / relative; before = path.read_bytes() if path.exists() else None
                mode = path.stat().st_mode if path.exists() else None
                if mode is not None: path.chmod(mode | 0o200)
                path.write_bytes(raw)
                with self.assertRaises(ValueError): B.check(self.args(imported=False))
                if before is None: path.unlink()
                else: path.write_bytes(before); path.chmod(mode)

    def test_missing_dependency_notice_refuses_ready_reuse_without_repairing_it(self):
        ready = B.install(self.args()); root = Path(ready["receipt"]).parent
        (root / "licenses/dependency/PANDAS_LICENSE").unlink()
        before = self.snapshot(self.cache)
        with self.assertRaises(FileNotFoundError): B.install(self.args())
        self.assertEqual(self.snapshot(self.cache), before)

    def test_changed_staged_notice_blocks_interrupted_publication(self):
        target = self.interrupt_publication()
        ledger, _, _ = B.load(B.publication_path(target))
        stage = target.parent / ledger['stage']
        (stage / 'licenses/EXPAT_LICENSE').write_bytes(b'changed prepared notice')
        before = self.snapshot(self.source)
        with self.assertRaises(ValueError): B.install(self.args())
        self.assertFalse((target / 'receipt.json').exists())
        self.assertEqual((stage / 'licenses/EXPAT_LICENSE').read_bytes(), b'changed prepared notice')
        self.assertEqual(self.snapshot(self.source), before)

    def test_different_explicit_evidence_cannot_silently_reuse_old_receipt(self):
        result = B.install(self.args()); before = Path(result["receipt"]).read_bytes()
        self.result["scope"] = "a different receipt"; self.refresh()
        with self.assertRaises(ValueError): B.install(self.args())
        self.assertEqual(Path(result["receipt"]).read_bytes(), before)

    def test_multiple_matching_builds_require_explicit_identity(self):
        first = B.install(self.args())
        self.binary.write_bytes(elf(b"second")); self.result["binary_sha256"] = B.digest(self.binary); self.refresh()
        B.install(self.args())
        with self.assertRaises(ValueError): B.check(self.args(imported=False))
        self.assertEqual(B.check(self.args(imported=False, identity=first["identity"])), first)

    def test_import_qualification_cannot_claim_different_producer_source(self):
        ready = B.install(self.args()); root = Path(ready["receipt"]).parent
        pins, _, _ = B.load(root / "evidence/build-pins.json")
        pins["source_files"]["rust/Cargo.lock"] = "0" * 64
        write_json(root / "evidence/build-pins.json", pins)
        result, _, _ = B.load(root / "evidence/validation.json")
        result["source_pins_sha256"] = B.digest(root / "evidence/build-pins.json")
        write_json(root / "evidence/validation.json", result)
        receipt, _, _ = B.load(root / "receipt.json")
        receipt["evidence_sha256"] = {name: B.digest(root / "evidence" / name) for name in receipt["evidence_sha256"]}
        write_json(root / "receipt.json", receipt)
        with self.assertRaises(ValueError): B.check(self.args(imported=False))

    def test_locked_cold_build_isolated_flags_and_no_inherited_qualification(self):
        state = B.source_state(self.source)
        tools = self.folder / "tools"; tools.mkdir()
        compiler = dict(channel="1.90.0", host=B.HOST, files={name: B.sha(name.encode()) for name in ("cargo", "rustc", "rustdoc", "rustfmt")})
        calls = []
        def cargo(argv, **kwargs):
            calls.append((argv, kwargs)); env = kwargs["env"]
            # include_str! in rust/src/cadd_hub.rs reads this file at compile time.
            self.assertEqual((kwargs["cwd"] / "../data/liftover/GRC_to_ucsc.contig.map.txt").read_text(), "1 chr1\n")
            path = Path(env["CARGO_TARGET_DIR"]) / "release/priva"
            path.parent.mkdir(parents=True); path.write_bytes(elf(b"cold")); path.chmod(0o700)
            kwargs["stdout"].write(b"synthetic compiler log\n")
        with mock.patch.object(B, "obtain_toolchain", return_value=tools), \
             mock.patch.object(B, "toolchain_descriptor", return_value=compiler), \
             mock.patch.object(B.subprocess, "run", side_effect=cargo), \
             mock.patch.dict(os.environ, {"RUSTFLAGS": "-C target-cpu=native", "RUSTC_WRAPPER": "foreign", "CARGO_PROFILE_RELEASE_OPT_LEVEL": "1"}):
            ready = B.install(self.args(imported=False, jobs=2))
        argv, kwargs = calls[0]
        self.assertIn("--locked", argv); self.assertIn("--release", argv)
        self.assertEqual(argv[argv.index("--jobs") + 1], "2")
        self.assertEqual(kwargs["env"]["RUSTFLAGS"], "")
        self.assertNotIn("RUSTC_WRAPPER", kwargs["env"])
        self.assertNotIn("CARGO_PROFILE_RELEASE_OPT_LEVEL", kwargs["env"])
        self.assertEqual(ready["qualification"], dict(origin="cold_locked_build", scientific_acceptance=False, production_backend_selected=False))
        self.assertEqual(B.source_state(self.source), state)
        self.assertEqual(B.check(self.args(imported=False)), ready)

    def test_bootstrap_download_checksum_is_checked_before_execution(self):
        path = self.folder / "rustup-init"
        with mock.patch.object(B.urllib.request, "urlopen", return_value=io.BytesIO(b"wrong")):
            with self.assertRaises(ValueError): B.download_bootstrap(path)
        self.assertEqual(path.read_bytes(), b"wrong")

    def test_xml_build_pins_static_bytes_and_preserves_explicit_search_path(self):
        prefix = self.folder / 'python'; pc = prefix / 'lib/pkgconfig'
        pc.mkdir(parents=True)
        archive = prefix / 'lib/libexpat.a'; archive.write_bytes(b'static library fixture')
        state = {'files': {'rust/src/interpro_mapping_builder/xml.rs': 'fixture'}}
        env = dict(PATH='/fixture/bin', PKG_CONFIG_PATH='/explicit/pkgconfig')
        with mock.patch.object(B.sys, 'prefix', str(prefix)), \
             mock.patch.object(B.shutil, 'which', return_value='/fixture/bin/pkg-config'), \
             mock.patch.object(B.subprocess, 'check_output', side_effect=['2.8.5\n', str(prefix/'lib')+'\n']):
            result = B.xml_build_dependency(state, env)
        self.assertEqual(result, {'name': 'expat', 'version': '2.8.5', 'static_archive_sha256': B.sha(archive.read_bytes())})
        self.assertEqual(env['PKG_CONFIG_PATH'], '/explicit/pkgconfig'+os.pathsep+str(pc))

    def test_xml_build_refuses_dynamic_only_environment(self):
        lib = self.folder / 'lib'; lib.mkdir(); (lib/'libexpat.so').write_bytes(b'dynamic fixture')
        state = {'files': {'rust/src/interpro_mapping_builder/xml.rs': 'fixture'}}
        with mock.patch.object(B.shutil, 'which', return_value='/fixture/pkg-config'), \
             mock.patch.object(B.subprocess, 'check_output', side_effect=['2.8.5\n', str(lib)+'\n']):
            with self.assertRaisesRegex(ValueError, 'static Expat archive'):
                B.xml_build_dependency(state, dict(PATH='/fixture'))

    def test_xml_dependency_change_during_compile_prevents_publication(self):
        state = B.source_state(self.source)
        tools = self.folder / 'tools'; tools.mkdir()
        def cargo(argv, **kwargs):
            binary = Path(kwargs['env']['CARGO_TARGET_DIR'])/'release/priva'
            binary.parent.mkdir(parents=True); binary.write_bytes(elf(b'cold')); binary.chmod(0o700)
        with mock.patch.object(B, 'obtain_toolchain', return_value=tools), \
             mock.patch.object(B, 'toolchain_descriptor', return_value={}), \
             mock.patch.object(B, 'xml_build_dependency', side_effect=[{'static_archive_sha256':'before'}, {'static_archive_sha256':'changed'}]), \
             mock.patch.object(B.subprocess, 'run', side_effect=cargo):
            with self.assertRaisesRegex(ValueError, 'dependency changed'):
                B.cold_build(self.cache, self.source, state, self.args(imported=False))
        self.assertFalse(list(self.cache.glob('*/receipt.json')))

    def test_automatic_toolchain_setup_uses_owned_roots_and_no_modify_path(self):
        state = B.source_state(self.source); base = B.managed_base(self.cache, create=True)
        compiler = dict(channel=state["channel"], host=B.HOST, files={})
        calls = []
        def run(argv, **kwargs):
            calls.append((argv, kwargs))
            if "--no-modify-path" in argv:
                stage = Path(kwargs["env"]["RUSTUP_HOME"]).parent
                (stage / "cargo/bin").mkdir(parents=True); (stage / "cargo/bin/rustup").write_bytes(b"fixture")
                (stage / "rustup/toolchains" / (state["channel"] + "-" + B.HOST) / "bin").mkdir(parents=True)
        with mock.patch.object(B.shutil, "which", return_value=None), \
             mock.patch.object(B.Path, "home", return_value=self.folder / "absent-home"), \
             mock.patch.dict(os.environ, {"RUSTUP_HOME": str(self.folder / "absent-rustup")}), \
             mock.patch.object(B, "download_bootstrap", side_effect=lambda path: path.write_bytes(b"fixture")), \
             mock.patch.object(B, "toolchain_descriptor", return_value=compiler), \
             mock.patch.object(B.subprocess, "run", side_effect=run):
            path = B.obtain_toolchain(base, state)
        self.assertTrue(path.is_relative_to(base))
        self.assertIn("--no-modify-path", calls[0][0])
        self.assertIn("--default-toolchain", calls[0][0])
        self.assertEqual(calls[0][0][-1], "1.90.0")
        self.assertTrue(Path(calls[0][1]["env"]["CARGO_HOME"]).is_relative_to(base))
        self.assertEqual(calls[1][0][1:4], ["component", "add", "rustfmt"])

    def test_incomplete_external_toolchain_is_preserved_and_falls_back(self):
        state = B.source_state(self.source); base = B.managed_base(self.cache, create=True)
        home = self.folder/'external-rustup'
        binary = home/'toolchains'/('1.90.0-'+B.HOST)/'bin/rustc'
        binary.parent.mkdir(parents=True); binary.write_bytes(b'incomplete external compiler')
        before = self.snapshot(home)
        with mock.patch.object(B.shutil, 'which', return_value=None), \
             mock.patch.object(B.Path, 'home', return_value=self.folder/'absent-home'), \
             mock.patch.dict(os.environ, {'RUSTUP_HOME':str(home)}), \
             mock.patch.object(B, 'toolchain_descriptor', side_effect=ValueError('incomplete external toolchain')):
            with self.assertRaisesRegex(ValueError, 'offline mode'):
                B.obtain_toolchain(base,state,offline=True)
        self.assertEqual(self.snapshot(home), before)

    def test_interrupted_owned_toolchain_setup_retries_without_ready_marker(self):
        state = B.source_state(self.source); base = B.managed_base(self.cache, create=True)
        compiler = dict(channel=state['channel'],host=B.HOST,files={})
        calls=[]
        def run(argv, **kwargs):
            calls.append(argv)
            target=Path(kwargs['env']['RUSTUP_HOME']).parent
            (target/'cargo/bin').mkdir(parents=True,exist_ok=True)
            (target/'cargo/bin/rustup').write_bytes(b'fixture')
            (target/'rustup/toolchains'/('1.90.0-'+B.HOST)/'bin').mkdir(parents=True,exist_ok=True)
            if len(calls)==1: raise subprocess.CalledProcessError(1,argv)
        with mock.patch.object(B.shutil,'which',return_value=None), \
             mock.patch.object(B.Path,'home',return_value=self.folder/'absent-home'), \
             mock.patch.dict(os.environ,{'RUSTUP_HOME':str(self.folder/'absent-rustup')}), \
             mock.patch.object(B,'RUSTUP_SHA256',B.sha(b'fixture')), \
             mock.patch.object(B,'download_bootstrap',side_effect=lambda path:path.write_bytes(b'fixture')) as download, \
             mock.patch.object(B,'toolchain_descriptor',return_value=compiler), \
             mock.patch.object(B.subprocess,'run',side_effect=run):
            with self.assertRaises(subprocess.CalledProcessError):B.obtain_toolchain(base,state)
            self.assertFalse((base/'toolchain-1.90.0/owner.json').exists())
            result=B.obtain_toolchain(base,state)
            self.assertTrue(result.is_dir());self.assertEqual(download.call_count,1)
        self.assertEqual(len(list((base/'toolchain-1.90.0').glob('install-*.log'))),2)


class ShellInstallTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory(dir=os.environ.get("PRIVA_TEST_TMPDIR", os.environ.get("TMPDIR")))
        self.addCleanup(self.temp.cleanup); self.folder = Path(self.temp.name)

    def run_shell(self, extra="", status=0, alias=False, binary=None):
        cache = self.folder / "cache"; cache.mkdir(exist_ok=True)
        selected = cache
        if alias:
            selected = self.folder / "alias"; selected.symlink_to(cache.name, target_is_directory=True)
        if binary is None:
            binary = self.folder / "managed-priva"; binary.write_text("#!/bin/sh\n"); binary.chmod(0o555)
        config = self.folder / "selected.yaml"
        config.write_text(f"cache_root: {selected}\npriva_binary: \"\"\n")
        calls = self.folder / "calls.json"
        script = r'''
source "$INSTALLER" >/dev/null || exit
read_yaml() { "$TEST_PYTHON" -c 'import sys,yaml;print(yaml.safe_load(open(sys.argv[1]))[sys.argv[2]])' "$1" "$2"; }
python() { "$TEST_PYTHON" -c 'import json,os,sys;open(os.environ["CALLS"],"w").write(json.dumps(sys.argv[1:]));print(json.dumps({"binary": os.environ["BINARY"], "ready": True}))' "$@"; return "$HELPER_STATUS"; }
PriVA_rust_install --config "$CONFIG" ''' + extra
        env = {key: value for key, value in os.environ.items() if not key.startswith("BASH_FUNC_") and key not in ("BASH_ENV", "ENV")}
        env.update(INSTALLER=str(ROOT / "scripts/install_utils.sh"), TEST_PYTHON=sys.executable,
                   CONFIG=str(config), CALLS=str(calls), HELPER_STATUS=str(status), BINARY=str(binary))
        result = subprocess.run(["bash", "-c", script], env=env, capture_output=True, text=True)
        return result, json.loads(calls.read_text()) if calls.exists() else None, cache, config

    def test_actual_shell_function_forwards_only_current_source_and_normalized_cache(self):
        result, calls, cache, config = self.run_shell(alias=True)
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertEqual(calls, [str(ROOT / "scripts/bootstrap_priva_rust.py"), "install", "--source-root", str(ROOT), "--cache-root", str(cache)])
        # The pipeline scripts find the program through this key.
        self.assertEqual(config.read_text().splitlines()[1], f"priva_binary: {self.folder / 'managed-priva'}")

    def test_helper_failure_propagates(self):
        result, _, _, config = self.run_shell(status=73)
        self.assertEqual(result.returncode, 1)
        self.assertIn('priva_binary: ""', config.read_text())

    def test_reported_path_must_be_an_executable(self):
        missing = self.folder / "missing-priva"
        result, _, _, config = self.run_shell(binary=missing)
        self.assertEqual(result.returncode, 1)
        self.assertIn('priva_binary: ""', config.read_text())

    def test_extra_arguments_refused_before_helper(self):
        result, calls, _, _ = self.run_shell(extra="--cache-root /foreign")
        self.assertEqual(result.returncode, 2); self.assertIsNone(calls)

    def test_every_compiled_in_file_is_copied_by_the_cold_build(self):
        # include_str!/include_bytes! outside the test modules must name files the
        # cold build copies (rust/src or SOURCE_FIXED); otherwise the build stops.
        import re
        listed = set(B.source_state(ROOT)["files"])
        for path in sorted((ROOT / "rust/src").rglob("*.rs")):
            code = path.read_text().split("#[cfg(test)]\nmod tests {", 1)[0]
            for name in re.findall(r'include_(?:str|bytes)!\("([^"]+)"\)', code):
                target = os.path.normpath(os.path.join(os.path.dirname(path.relative_to(ROOT)), name))
                self.assertIn(target, listed, f"{path.relative_to(ROOT)} compiles in {name}")

    def test_main_calls_rust_once_after_driver_without_backend_mutation(self):
        source = (ROOT / "scripts/install_utils.sh").read_text()
        # The entry point's subshell contains its private builder-view marker
        # and lock descriptor, including early returns.
        main = source.split("function main_install() (", 1)[1].split("\n)\n", 1)[0]
        self.assertEqual(main.count('PriVA_rust_install --config "${config_file}"'), 1)
        self.assertLess(main.index('Nextflow_runtime_install "${config_file}"'), main.index('PriVA_rust_install --config'))
        self.assertLess(main.index('PriVA_rust_install --config'), main.index('reference_genome_install "${config_file}"'))
        function = source.split("function PriVA_rust_install() {", 1)[1].split("\n}\n", 1)[0]
        # It saves the program's path as priva_binary and nothing else: no ranking
        # selection, environment or PATH change.
        self.assertEqual(function.count('update_or_append_yaml'), 1)
        self.assertIn('update_or_append_yaml "$2" priva_binary "${binary}"', function)
        for forbidden in ('priva_ranking_binary', 'PRIVA_RUST', 'export ', 'PATH='):
            self.assertNotIn(forbidden, function)


if __name__ == "__main__":
    unittest.main()
