"""Owned VEP/assembly recovery, exact retention, and retry barriers on tiny files.

Scientific generation is injected; existing installer/candidate fixtures exercise
the unchanged stage argv. No production cache or network resource is accessed.
"""
from __future__ import annotations

import gzip
import json
import os
from pathlib import Path
import re
import subprocess
import sys
import tempfile
import unittest
from unittest import mock

import yaml

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
import alphamissense_cache as alpha
import clinvar_gene_refresh as refresh


def write_vcf(path, *, major="116", assembly="GRCh38", headers=None, label="old"):
    path.parent.mkdir(parents=True, exist_ok=True)
    if headers is None:
        headers = [f'##VEP="v{major}.3" assembly="{assembly}"']
    body = ("##fileformat=VCFv4.2\n" + "\n".join(headers) + "\n"
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
            f"chr1\t1\t{label}\tA\tC\t.\tPASS\t.\n")
    path.write_bytes(gzip.compress(body.encode()) if path.name.endswith(".gz") else body.encode())
    return path


class ClinVarOwnedCacheTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory(prefix="clinvar_owned_", dir=os.environ.get("PRIVA_TEST_TMPDIR"))
        self.addCleanup(self.temp.cleanup)
        self.folder = Path(self.temp.name)
        self.cache = self.folder / "cache"
        self.vcf_dir = self.cache / "clinvar"
        self.target = self.vcf_dir / "clinvar.hg38.vep.vcf.gz"
        self.config = self.folder / "config.yaml"
        self.settings = dict(cache_root=str(self.cache), clinvar_vcf_dir=str(self.vcf_dir),
                             clinvar_vcf=str(self.target), assembly="hg38", vep_installed_version="116")
        self.save_config()
        self.calls = []

    def save_config(self):
        self.config.write_text(yaml.safe_dump(self.settings))

    def context(self):
        return refresh.deployment_context(self.config)

    def rebuild(self, *arguments):
        self.calls.append(arguments)
        write_vcf(self.target, label="new")
        Path(str(self.target) + ".tbi").write_bytes(b"new index")

    def deploy(self, runner=None, **kwargs):
        return refresh.deploy_cache(self.config, ROOT, runner=self.rebuild if runner is None else runner, **kwargs)

    def journals(self):
        return [json.loads(path.read_text()) for path in self.cache.rglob("retirement.json")]

    def assert_retained(self, original, content):
        moves = [item for journal in self.journals() for item in journal["moves"]
                 if item["source"] == str(original) and item["moved"]]
        self.assertTrue(moves, str(original))
        self.assertTrue(any(Path(item["retained"]).read_bytes() == content for item in moves))

    def test_matching_headers_reuse_without_retirement_and_preserve_source(self):
        for extra in ("", "##INFO=<ID=AF_joint,Number=A,Type=Float,Description=\"fixture\">"):
            write_vcf(self.target, headers=['##VEP="v116.3" assembly="GRCh38.p14"'] + ([extra] if extra else []))
            before = refresh.retirement_identity(self.target), self.target.read_bytes()
            seen = []
            self.deploy(runner=lambda *args: seen.append(args))
            self.assertEqual(len(seen), 1)
            self.assertEqual(seen[0][1], "_ClinVar_VCF_deploy")
            self.assertEqual((refresh.retirement_identity(self.target), self.target.read_bytes()), before)
            self.assertEqual(self.journals(), [])

    def test_owned_old_wrong_missing_and_ambiguous_metadata_retain_then_rebuild(self):
        for headers in (['##VEP="v113" assembly="GRCh38"'], ['##VEP="v116" assembly="GRCh37"'],
                        [], ['##VEP="v116" assembly="GRCh38"'] * 2):
            with self.subTest(headers=headers):
                write_vcf(self.target, headers=headers)
                index = Path(str(self.target) + ".tbi")
                index.write_bytes(b"old index")
                ready = refresh.candidate_receipt_path(self.target)
                ready.write_text('{"status":"old evidence"}\n')
                before = {path: path.read_bytes() for path in (self.target, index, ready)}
                unrelated = self.vcf_dir / "unrelated.keep"
                unrelated.write_bytes(b"unrelated retained in place")
                def rebuild(*args):
                    self.assertTrue(all(not path.exists() for path in before))
                    self.rebuild(*args)
                self.deploy(runner=rebuild)
                for path, content in before.items():
                    self.assert_retained(path, content)
                self.assertEqual(unrelated.read_bytes(), b"unrelated retained in place")
                self.assertFalse(refresh.deployment_pending_path(self.context()).exists())
        self.assertTrue(all(value["scientific_qualified"] is False for value in self.journals()))
        self.assertTrue(all("old_profile" not in value for value in self.journals()))

    def test_explicit_incompatible_or_unknown_input_is_read_only(self):
        custom = self.folder / "supplied.vcf.gz"
        self.settings["clinvar_vcf"] = str(custom)
        self.save_config()
        for headers in ([], ['##VEP="v113" assembly="GRCh38"'], ['##VEP="v116" assembly="GRCh37"']):
            write_vcf(custom, headers=headers)
            before = custom.read_bytes(), refresh.retirement_identity(custom), self.config.read_bytes()
            with self.assertRaisesRegex(ValueError, "Explicit ClinVar annotation.*preserved"):
                self.deploy()
            self.assertEqual((custom.read_bytes(), refresh.retirement_identity(custom), self.config.read_bytes()), before)
        self.assertEqual(self.calls, [])
        self.assertFalse(self.cache.exists())

    def test_explicit_compatible_input_is_not_redirected_or_modified(self):
        custom = write_vcf(self.folder / "supplied.vcf.gz")
        self.settings["clinvar_vcf"] = str(custom)
        self.save_config()
        before = custom.read_bytes(), self.config.read_bytes()
        seen = []
        self.deploy(runner=lambda *args: seen.append(args))
        self.assertEqual(len(seen), 1)
        self.assertEqual((custom.read_bytes(), self.config.read_bytes()), before)
        self.assertFalse(self.cache.exists())

    def test_linked_target_parent_escape_and_linked_companion_fail_before_moves(self):
        external = write_vcf(self.folder / "external.vcf.gz", major="113")
        self.vcf_dir.mkdir(parents=True)
        self.target.symlink_to(external)
        with self.assertRaisesRegex(ValueError, "Linked ClinVar"):
            self.deploy()
        self.target.unlink()
        self.vcf_dir.rmdir()
        outside = self.folder / "outside"
        outside.mkdir()
        self.vcf_dir.symlink_to(outside, target_is_directory=True)
        with self.assertRaisesRegex(ValueError, "inside.*cache_root"):
            self.deploy()
        self.vcf_dir.unlink()
        write_vcf(self.target, major="113")
        index = Path(str(self.target) + ".tbi")
        index.symlink_to(external)
        with self.assertRaisesRegex(ValueError, "Nonregular/linked"):
            self.deploy()
        self.assertEqual(self.calls, [])
        self.assertEqual(self.journals(), [])
        self.assertEqual(external.read_bytes(), self.target.read_bytes())

    def test_config_and_pinned_source_aliases_fail_without_retention_or_lock(self):
        write_vcf(self.target, major="113")
        raw = self.vcf_dir / "clinvar.GRCh38.vcf.gz"
        raw.write_bytes(b"pinned raw")
        before = self.target.read_bytes(), raw.read_bytes(), self.config.read_bytes()
        with self.assertRaisesRegex(ValueError, "selected config/source"):
            self.deploy(source_vcf=raw)
        self.assertEqual((self.target.read_bytes(), raw.read_bytes(), self.config.read_bytes()), before)
        self.assertFalse((self.vcf_dir / ".clinvar-deploy-locks").exists())
        raw.unlink()
        index = Path(str(self.target) + ".tbi")
        os.link(self.config, index)
        with self.assertRaisesRegex(ValueError, "selected config/source"):
            self.deploy()
        self.assertEqual(self.journals(), [])

    def test_pinned_noncanonical_shard_namespace_and_hardlink_alias_fail_before_retention(self):
        write_vcf(self.target, major="113")
        shard = self.vcf_dir / "clinvar.hg38.vep.chr1_KI_random.vcf.gz"
        shard.write_bytes(b"pinned source in shard namespace")
        before = shard.read_bytes(), self.target.read_bytes()
        with self.assertRaisesRegex(ValueError, "may alias a managed overwrite"):
            self.deploy(source_vcf=shard)
        separate = self.folder / "sources/raw.vcf.gz"
        separate.parent.mkdir()
        os.link(shard, separate)
        with self.assertRaisesRegex(ValueError, "may alias a managed overwrite"):
            self.deploy(source_vcf=separate)
        self.assertEqual((shard.read_bytes(), self.target.read_bytes()), before)
        self.assertEqual(self.journals(), [])
        self.assertFalse((self.vcf_dir / ".clinvar-deploy-locks").exists())

    def test_separate_single_link_pinned_source_and_index_are_preserved(self):
        write_vcf(self.target, major="113")
        source = write_vcf(self.folder / "sources/raw.vcf.gz")
        source_index = Path(str(source) + ".tbi")
        source_index.write_bytes(b"pinned index")
        before = [(path.read_bytes(), refresh.retirement_identity(path)) for path in (source, source_index)]
        def rebuild(*args):
            self.assertEqual(list(args[-2:]), ["--source-vcf", source])
            self.rebuild(*args)
            refresh.retire_rebuilt_shards(self.config, self.target, ["chr1"])
        self.deploy(runner=rebuild, source_vcf=source)
        self.assertEqual([(path.read_bytes(), refresh.retirement_identity(path)) for path in (source, source_index)], before)

    def test_changed_header_never_authorizes_retirement_even_on_parse_or_mismatch(self):
        original = alpha._check_vep_annotation_header
        for headers in (None, [], ['##VEP="v113" assembly="GRCh38"']):
            write_vcf(self.target, headers=headers)
            def change(*args):
                try:
                    return original(*args)
                finally:
                    self.target.write_bytes(self.target.read_bytes() + b"mutation")
            with mock.patch.object(alpha, "_check_vep_annotation_header", side_effect=change):
                with self.assertRaises(alpha.AnnotationHeaderChangedError):
                    self.deploy()
            self.assertEqual(self.journals(), [])
            self.assertEqual(self.calls, [])

    def test_change_between_failed_header_and_retirement_is_not_retired(self):
        write_vcf(self.target, major="113")
        original = refresh.retire_owned_cache
        def changed(*args, **kwargs):
            write_vcf(self.target, major="117", label="concurrent")
            return original(*args, **kwargs)
        with mock.patch.object(refresh, "retire_owned_cache", side_effect=changed):
            with self.assertRaises(alpha.AnnotationHeaderChangedError):
                self.deploy()
        self.assertEqual(self.journals(), [])
        self.assertEqual(self.calls, [])

    def test_newer_wrong_profile_plain_repair_is_retained_before_body_can_replace_gzip(self):
        write_vcf(self.target)
        plain = write_vcf(Path(str(self.target)[:-3]), major="113")
        old = {path: path.read_bytes() for path in (self.target, plain)}
        def rebuild(*args):
            self.assertFalse(plain.exists())
            self.assertFalse(self.target.exists())
            self.rebuild(*args)
        self.deploy(runner=rebuild)
        for path, content in old.items():
            self.assert_retained(path, content)

    def test_primary_partial_move_failure_is_retained_and_retry_rebuilds(self):
        write_vcf(self.target, major="113")
        index = Path(str(self.target) + ".tbi")
        index.write_bytes(b"index")
        ready = refresh.candidate_receipt_path(self.target)
        ready.write_bytes(b"ready evidence")
        old = {path: path.read_bytes() for path in (self.target, index, ready)}
        rename = Path.rename
        def fail_second(path, output):
            if path == ready:
                raise OSError("injected second move failure")
            return rename(path, output)
        with mock.patch.object(Path, "rename", new=fail_second):
            with self.assertRaisesRegex(OSError, "second move"):
                self.deploy()
        self.assertEqual(self.calls, [])
        self.assertEqual(self.target.read_bytes(), old[self.target])
        self.assertEqual(self.journals()[0]["status"], "retirement_incomplete")
        self.deploy()
        for path, content in old.items():
            self.assert_retained(path, content)
        self.assertEqual(len(self.calls), 1)

    def test_failed_generator_keeps_old_journal_and_new_partial_for_retry(self):
        write_vcf(self.target, major="113")
        original = self.target.read_bytes()
        def failed(*args):
            self.target.write_bytes(b"partial VEP output")
            raise RuntimeError("injected VEP failure")
        with self.assertRaisesRegex(RuntimeError, "VEP failure"):
            self.deploy(runner=failed)
        self.assertEqual(self.target.read_bytes(), b"partial VEP output")
        self.assert_retained(self.target, original)
        self.deploy()
        self.assert_retained(self.target, b"partial VEP output")

    def test_journal_failure_before_first_move_never_starts_rebuild_and_retry_retains(self):
        write_vcf(self.target, major="113")
        before = self.target.read_bytes()
        original = refresh.write_json
        def fail_journal(path, value):
            if path.name == "retirement.json":
                raise OSError("injected retirement journal failure")
            return original(path, value)
        with mock.patch.object(refresh, "write_json", side_effect=fail_journal):
            with self.assertRaisesRegex(OSError, "journal failure"):
                self.deploy()
        self.assertEqual(self.target.read_bytes(), before)
        self.assertEqual(self.calls, [])
        self.deploy()
        self.assert_retained(self.target, before)

    def setup_shard_recovery(self, *, source_vcf=None):
        write_vcf(self.target, major="113")
        refresh.retire_owned_cache(self.config, reason="test old generation", source_vcf=source_vcf)
        write_vcf(self.target, label="new VEP")
        names = ["chr1", "chr1_KI270706v1_random", "prefixchr1suffix"]
        paths = refresh.shard_retirement_members(self.context(), names)
        for path in paths:
            path.write_bytes((path.name + " old bytes").encode())
            os.utime(path, ns=(2100000000000000000, 2100000000000000000))
        return names, {path: path.read_bytes() for path in paths}

    def pinned_source_fixture(self):
        source = write_vcf(self.folder / "sources/raw.vcf.gz")
        paths = [source, Path(str(source) + ".tbi"), Path(str(source) + ".csi")]
        for path in paths[1:]:
            path.write_bytes((path.name + " pinned index").encode())
        return source, paths

    def other_native_device(self):
        original = refresh.retirement_identity
        def observed(path):
            value = original(path)
            value["dev"] = str(int(value["dev"]) + 1)
            return value
        return mock.patch.object(refresh, "retirement_identity", side_effect=observed)

    def assert_namespace_unchanged(self, before):
        self.assertEqual({path: path.read_bytes() for path in before}, before)

    def test_primary_partial_retry_accepts_only_device_change_and_retains_raw_observations(self):
        source, sources = self.pinned_source_fixture()
        write_vcf(self.target, major="113")
        Path(str(self.target) + ".tbi").write_bytes(b"old index")
        Path(str(self.target) + ".csi").write_bytes(b"old alternative index")
        rename = Path.rename
        moved = 0
        def interrupt(path, output):
            nonlocal moved
            moved += 1
            if moved == 2:
                raise OSError("injected primary interruption")
            return rename(path, output)
        with mock.patch.object(Path, "rename", new=interrupt):
            with self.assertRaisesRegex(OSError, "primary interruption"):
                refresh.retire_owned_cache(self.config, reason="old generation", source_vcf=source)
        historical = refresh.deployment_pending(self.context())["protected_sources"]
        prior_journals = {path: path.read_bytes() for path in self.cache.rglob("retirement.json")}
        source_before = {path: (path.read_bytes(), refresh.retirement_identity(path)) for path in sources}
        with self.other_native_device():
            refresh.retire_owned_cache(self.config, reason="resume old generation", source_vcf=source)
        pending = refresh.deployment_pending(self.context())
        self.assertEqual(pending["protected_sources"], historical)
        current = pending["retry_observations"][-1]["protected_sources"]
        for path, previous in historical.items():
            self.assertEqual(current[path], {**previous, "dev": str(int(previous["dev"]) + 1)})
        self.assertFalse(self.target.exists())
        self.assertEqual({path: (path.read_bytes(), refresh.retirement_identity(path)) for path in sources}, source_before)
        self.assert_namespace_unchanged(prior_journals)
        self.assertEqual(self.calls, [])

    def test_partial_shard_retry_accepts_new_device_without_rewriting_prior_target_or_sources(self):
        source, sources = self.pinned_source_fixture()
        names, old = self.setup_shard_recovery(source_vcf=source)
        target_before = self.target.read_bytes(), refresh.retirement_identity(self.target)
        rename = Path.rename
        failing = list(old)[3]
        def interrupt(path, output):
            if path == failing:
                raise OSError("injected cross-node shard interruption")
            return rename(path, output)
        with mock.patch.object(Path, "rename", new=interrupt):
            with self.assertRaisesRegex(OSError, "shard interruption"):
                refresh.retire_rebuilt_shards(self.config, self.target, names)
        historical = refresh.deployment_pending(self.context())
        self.assertEqual(historical["phase"], "retiring_shards")
        prior_journals = {path: path.read_bytes() for path in self.cache.rglob("retirement.json")}
        with self.other_native_device():
            refresh.retire_rebuilt_shards(self.config, self.target, names)
        pending = refresh.deployment_pending(self.context())
        self.assertEqual(pending["phase"], "shards_retained")
        self.assertEqual(pending["rebuilt_target"], historical["rebuilt_target"])
        self.assertEqual(pending["protected_sources"], historical["protected_sources"])
        observation = pending["retry_observations"][-1]
        self.assertEqual(observation["rebuilt_target"],
                         {**historical["rebuilt_target"], "dev": str(int(historical["rebuilt_target"]["dev"]) + 1)})
        for path, previous in historical["protected_sources"].items():
            self.assertEqual(observation["protected_sources"][path],
                             {**previous, "dev": str(int(previous["dev"]) + 1)})
        self.assertEqual((self.target.read_bytes(), refresh.retirement_identity(self.target)), target_before)
        for path, content in old.items():
            self.assertFalse(path.exists())
            self.assert_retained(path, content)
        self.assertTrue(all(path.exists() for path in sources))
        self.assert_namespace_unchanged(prior_journals)
        self.assertEqual(self.calls, [])

    def test_completed_shard_gate_accepts_new_device_and_preserves_new_splits(self):
        source, _ = self.pinned_source_fixture()
        names, _ = self.setup_shard_recovery(source_vcf=source)
        refresh.retire_rebuilt_shards(self.config, self.target, names)
        historical = refresh.deployment_pending(self.context())
        new_shard = refresh.shard_retirement_members(self.context(), names)[0]
        new_shard.write_bytes(b"completed current-generation split")
        before = {path: path.read_bytes() for path in (self.target, new_shard)}
        journals = len(self.journals())
        with self.other_native_device():
            refresh.retire_rebuilt_shards(self.config, self.target, names)
        pending = refresh.deployment_pending(self.context())
        self.assertEqual(pending["rebuilt_target"], historical["rebuilt_target"])
        self.assertEqual(pending["protected_sources"], historical["protected_sources"])
        self.assertEqual(len(pending["retry_observations"]), len(historical["retry_observations"]) + 1)
        self.assert_namespace_unchanged(before)
        self.assertEqual(len(self.journals()), journals)

    def test_persisted_target_requires_all_fields_and_rejects_every_other_change(self):
        names, old = self.setup_shard_recovery()
        with mock.patch.object(Path, "rename", side_effect=OSError("stop before first shard move")):
            with self.assertRaisesRegex(OSError, "first shard move"):
                refresh.retire_rebuilt_shards(self.config, self.target, names)
        fields = ("dev", "ino", "mode", "nlink", "size", "mtime_ns", "ctime_ns")
        mutations = [(field, "missing") for field in fields]
        mutations += [(field, "changed") for field in fields if field != "dev"]
        mutations += [("dev", "malformed"), ("extra", "unexpected")]
        for phase in ("retiring_shards", "shards_retained"):
            if phase == "shards_retained":
                refresh.retire_rebuilt_shards(self.config, self.target, names)
                list(old)[0].write_bytes(b"new-generation split must survive rejection")
            baseline = refresh.deployment_pending(self.context())
            self.assertEqual(baseline["phase"], phase)
            for field, operation in mutations:
                with self.subTest(phase=phase, field=field, operation=operation):
                    pending = json.loads(json.dumps(baseline))
                    recorded = pending["rebuilt_target"]
                    if operation == "missing":
                        recorded.pop(field)
                    elif operation == "changed":
                        recorded[field] = str(int(recorded[field]) + 1)
                    else:
                        recorded[field] = "not a native identity"
                    pending_path = refresh.deployment_pending_path(self.context())
                    refresh.write_json(pending_path, pending)
                    before = {path: path.read_bytes() for path in [self.target, pending_path, *old] if path.exists()}
                    with self.assertRaisesRegex(ValueError, "output/selected contigs changed"):
                        refresh.retire_rebuilt_shards(self.config, self.target, names)
                    self.assert_namespace_unchanged(before)
            refresh.write_json(refresh.deployment_pending_path(self.context()), baseline)

    def test_persisted_pinned_sources_require_all_fields_in_each_retry_phase(self):
        source, sources = self.pinned_source_fixture()
        write_vcf(self.target, major="113")
        refresh.retire_owned_cache(self.config, reason="old generation", source_vcf=source)
        fields = ("dev", "ino", "mode", "nlink", "size", "mtime_ns", "ctime_ns")
        mutations = [(field, "missing") for field in fields]
        mutations += [(field, "changed") for field in fields if field != "dev"]
        mutations += [("dev", "malformed"), ("extra", "unexpected")]
        for phase in ("awaiting_vep", "retiring_shards", "shards_retained"):
            if phase == "retiring_shards":
                write_vcf(self.target, label="new VEP")
                names = ["chr1", "chr1_KI270706v1_random"]
                for path in refresh.shard_retirement_members(self.context(), names):
                    path.write_bytes(b"old shard")
                with mock.patch.object(Path, "rename", side_effect=OSError("stop before first shard move")):
                    with self.assertRaisesRegex(OSError, "first shard move"):
                        refresh.retire_rebuilt_shards(self.config, self.target, names)
            elif phase == "shards_retained":
                refresh.retire_rebuilt_shards(self.config, self.target, names)
            baseline = refresh.deployment_pending(self.context())
            self.assertEqual(baseline["phase"], phase)
            for path in sources:
                for field, operation in mutations:
                    with self.subTest(phase=phase, source=path.name, field=field, operation=operation):
                        pending = json.loads(json.dumps(baseline))
                        recorded = pending["protected_sources"][str(path)]
                        if operation == "missing":
                            recorded.pop(field)
                        elif operation == "changed":
                            recorded[field] = str(int(recorded[field]) + 1)
                        else:
                            recorded[field] = "not a native identity"
                        pending_path = refresh.deployment_pending_path(self.context())
                        refresh.write_json(pending_path, pending)
                        before = {item: item.read_bytes() for item in [pending_path, *sources, self.target] if item.exists()}
                        with self.assertRaisesRegex(ValueError, "Pinned ClinVar source/index.*changed"):
                            if phase == "awaiting_vep":
                                refresh.retire_owned_cache(self.config, reason="retry", source_vcf=source)
                            else:
                                refresh.retire_rebuilt_shards(self.config, self.target, names)
                        self.assert_namespace_unchanged(before)
            refresh.write_json(refresh.deployment_pending_path(self.context()), baseline)

    def test_primary_retry_keeps_exact_pinned_source_path_set(self):
        source, sources = self.pinned_source_fixture()
        write_vcf(self.target, major="113")
        refresh.retire_owned_cache(self.config, reason="old generation", source_vcf=source)
        baseline = refresh.deployment_pending(self.context())
        for operation in ("missing", "extra"):
            with self.subTest(operation=operation):
                pending = json.loads(json.dumps(baseline))
                if operation == "missing":
                    pending["protected_sources"].pop(str(sources[-1]))
                else:
                    pending["protected_sources"][str(source) + ".other"] = pending["protected_sources"][str(source)]
                pending_path = refresh.deployment_pending_path(self.context())
                refresh.write_json(pending_path, pending)
                before = {path: path.read_bytes() for path in [pending_path, *sources]}
                with self.assertRaisesRegex(ValueError, "identity changed during recovery"):
                    with self.other_native_device():
                        refresh.retire_owned_cache(self.config, reason="retry", source_vcf=source)
                self.assert_namespace_unchanged(before)

    def test_same_process_device_changes_still_stop_header_and_move_checks(self):
        names, _ = self.setup_shard_recovery()
        original = refresh.retirement_identity
        count = 0
        def during_header(path):
            nonlocal count
            value = original(path)
            if path == self.target:
                count += 1
                if count > 1:
                    value["dev"] = str(int(value["dev"]) + 1)
            return value
        before = self.target.read_bytes()
        with mock.patch.object(refresh, "retirement_identity", side_effect=during_header):
            with self.assertRaisesRegex(ValueError, "changed before shard retirement"):
                refresh.retire_rebuilt_shards(self.config, self.target, names)
        self.assertEqual(self.target.read_bytes(), before)
        def during_move(path):
            value = original(path)
            if ".installer-retained" in path.parts:
                value["dev"] = str(int(value["dev"]) + 1)
            return value
        with mock.patch.object(refresh, "retirement_identity", side_effect=during_move):
            with self.assertRaisesRegex(ValueError, "changed while being retained"):
                refresh.retire_rebuilt_shards(self.config, self.target, names)
        self.assertEqual(self.target.read_bytes(), before)
        self.assertEqual(self.calls, [])

    def test_interruption_after_vep_preserves_new_output_and_retires_all_literal_shards(self):
        names, old = self.setup_shard_recovery()
        new_before = self.target.read_bytes(), refresh.retirement_identity(self.target)
        unrelated = self.vcf_dir / "clinvar.hg38.vep.chr2.vcf.gz"
        unrelated.write_bytes(b"unselected untouched")
        refresh.retire_rebuilt_shards(self.config, self.target, names)
        self.assertEqual((self.target.read_bytes(), refresh.retirement_identity(self.target)), new_before)
        self.assertEqual(refresh.deployment_pending(self.context())["phase"], "shards_retained")
        for path, content in old.items():
            self.assertFalse(path.exists(), str(path))
            self.assert_retained(path, content)
        self.assertEqual(unrelated.read_bytes(), b"unselected untouched")
        self.assertEqual(len(self.journals()), 2)

    def test_partial_shard_retirement_retry_never_retires_new_vep(self):
        names, old = self.setup_shard_recovery()
        new_before = self.target.read_bytes(), refresh.retirement_identity(self.target)
        failing = list(old)[3]
        rename = Path.rename
        def failed(path, output):
            if path == failing:
                raise OSError("injected shard move failure")
            return rename(path, output)
        with mock.patch.object(Path, "rename", new=failed):
            with self.assertRaisesRegex(OSError, "shard move"):
                refresh.retire_rebuilt_shards(self.config, self.target, names)
        self.assertEqual(refresh.deployment_pending(self.context())["phase"], "retiring_shards")
        with self.assertRaisesRegex(ValueError, "new VEP output was preserved"):
            refresh.retire_owned_cache(self.config, reason="must not restart primary")
        refresh.retire_rebuilt_shards(self.config, self.target, names)
        self.assertEqual((self.target.read_bytes(), refresh.retirement_identity(self.target)), new_before)
        for path, content in old.items():
            self.assertFalse(path.exists())
            self.assert_retained(path, content)
        self.assertEqual(sum(journal["phase"] == "primary" for journal in self.journals()), 1)

    def test_context_changed_target_or_contigs_never_reuse_prior_shard_authorization(self):
        names, old = self.setup_shard_recovery()
        wrong = self.vcf_dir / "other.vcf.gz"
        with self.assertRaisesRegex(ValueError, "context does not match"):
            refresh.retire_rebuilt_shards(self.config, wrong, names)
        self.assertTrue(all(path.read_bytes() == content for path, content in old.items()))
        refresh.retire_rebuilt_shards(self.config, self.target, names)
        with self.assertRaisesRegex(ValueError, "output/selected contigs changed"):
            refresh.retire_rebuilt_shards(self.config, self.target, names[:-1])
        write_vcf(self.target, label="changed after retention")
        with self.assertRaisesRegex(ValueError, "output/selected contigs changed"):
            refresh.retire_rebuilt_shards(self.config, self.target, names)

    def test_path_escape_or_linked_shard_aborts_before_any_shard_move(self):
        names, old = self.setup_shard_recovery()
        with self.assertRaisesRegex(ValueError, "Unsafe literal"):
            refresh.retire_rebuilt_shards(self.config, self.target, names + ["chr1/../../escape"])
        linked = list(old)[-1]
        linked.unlink()
        linked.symlink_to(self.config)
        with self.assertRaisesRegex(ValueError, "Nonregular/linked"):
            refresh.retire_rebuilt_shards(self.config, self.target, names)
        for path, content in old.items():
            if path != linked:
                self.assertEqual(path.read_bytes(), content)

    def run_gnomad_until_split(self, *, context=None, input_vcf=None):
        text = (ROOT / "scripts/annotation_vcf.sh").read_text()
        function = re.search(r"^function anno_agg_gnomAD_data \(\) \{.*?^\}", text, re.M | re.S)[0]
        prelude = r'''
log() { :; }; randomID() { printf fixture; }
anno_agg_gnomAD_index() { return 0; }
check_vcf_validity() { return 0; }; check_vcf_infotags() { return 1; }
bcftools() { [[ $1 == query ]] || return 90; printf '%s\n' chr1 chr1_KI270706v1_random prefixchr1suffix chrM; }
parallel() { printf 'parallel-start\n' >&2; return 77; }
'''
        environment = dict(os.environ, SCRIPT_DIR=str(ROOT / "scripts"))
        environment.pop("PRIVA_CLINVAR_RECOVERY_CONFIG", None)
        if context is not None:
            environment["PRIVA_CLINVAR_RECOVERY_CONFIG"] = str(context)
        return subprocess.run(["bash", "-c", prelude + function + '\nanno_agg_gnomAD_data "$1" 1 hg38 unused',
                               "fixture", str(self.target if input_vcf is None else input_vcf)],
                              env=environment, capture_output=True, text=True, timeout=15)

    def test_no_context_original_hook_is_inert_and_context_passes_actual_literal_selection(self):
        names, old = self.setup_shard_recovery()
        result = self.run_gnomad_until_split()
        self.assertEqual(result.returncode, 77, result.stderr)
        self.assertIn("parallel-start", result.stderr)
        self.assertTrue(all(path.read_bytes() == content for path, content in old.items()))
        result = self.run_gnomad_until_split(context=self.config)
        self.assertEqual(result.returncode, 77, result.stderr)
        self.assertEqual(refresh.deployment_pending(self.context())["chromosomes"], names)
        self.assertTrue(all(not path.exists() for path in old))

    def test_hook_context_mismatch_fails_before_parallel_or_moves(self):
        _, old = self.setup_shard_recovery()
        result = self.run_gnomad_until_split(context=self.config, input_vcf=self.vcf_dir / "other.vcf.gz")
        self.assertNotEqual(result.returncode, 0)
        self.assertNotIn("parallel-start", result.stderr)
        self.assertIn("context does not match", result.stderr)
        self.assertTrue(all(path.read_bytes() == content for path, content in old.items()))

    def test_nested_refresh_workspace_lock_does_not_deadlock(self):
        write_vcf(self.target)
        code = """import sys; from pathlib import Path
sys.path.insert(0, sys.argv[1]); import clinvar_gene_refresh as r
config=Path(sys.argv[2]); folder=Path(sys.argv[3])
with r.build_lock(folder):
    r.deploy_cache(config, Path(sys.argv[4]), runner=lambda *args: None)
"""
        result = subprocess.run([sys.executable, "-c", code, str(ROOT / "scripts"), str(self.config),
                                 str(self.vcf_dir), str(ROOT)], capture_output=True, text=True, timeout=10)
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertTrue((self.vcf_dir / ".refresh.lock").exists())
        self.assertTrue((refresh.deployment_pending_path(self.context()).parent / ".refresh.lock").exists())

    def test_public_shell_forwards_config_and_named_source_without_changing_argv(self):
        text = (ROOT / "scripts/install_utils.sh").read_text()
        wrapper = re.search(r"^function ClinVar_VCF_deploy\(\) \{.*?^\}", text, re.M | re.S)[0]
        body = 'python() { printf "%s\\n" "$@"; };\n' + wrapper + '\nClinVar_VCF_deploy "$1" --source-vcf "$2"'
        source = self.folder / "input source.vcf.gz"
        result = subprocess.run(["bash", "-c", body, "fixture", str(self.config), str(source)],
            env=dict(os.environ, SCRIPT_DIR=str(ROOT / "scripts"), BASE_DIR=str(ROOT)), capture_output=True, text=True)
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertEqual(result.stdout.splitlines(), [str(ROOT / "scripts/clinvar_gene_refresh.py"), "deploy",
            "--config", str(self.config), "--source-root", str(ROOT), "--source-vcf", str(source)])


if __name__ == "__main__":
    unittest.main()
