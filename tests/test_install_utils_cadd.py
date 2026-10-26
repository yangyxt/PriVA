"""Tiny CADD cache fixtures only; no publisher downloads or real cache access."""
import contextlib
import errno
import hashlib
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer
import importlib.util
import io
import json
import os
from pathlib import Path
import subprocess
import sys
import tarfile
import tempfile
import threading
from types import SimpleNamespace
import unittest
from unittest import mock


ROOT = Path(__file__).resolve().parents[1]
SPEC = importlib.util.spec_from_file_location("cadd_cache_install", ROOT / "scripts/cadd_cache_install.py")
CADD = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(CADD)


def md5(value):
    return hashlib.md5(value).hexdigest()


class CaddCacheTests(unittest.TestCase):
    def setUp(self):
        temporary = tempfile.TemporaryDirectory(prefix="priva_cadd_", dir=os.environ.get("PRIVA_TEST_TMPDIR", os.environ.get("TMPDIR")))
        self.addCleanup(temporary.cleanup)
        self.directory = Path(temporary.name)
        self.sources = self.directory / "sources"
        self.sources.mkdir()
        self.base = self.directory / "cache/CADD-scripts-1.7.2"
        self.base.mkdir(parents=True)
        (self.base / "CADD.sh").write_text("#!/bin/bash\nexit 99\n")
        self.output = contextlib.redirect_stdout(io.StringIO())
        self.output.__enter__()
        self.addCleanup(self.output.__exit__, None, None, None)
        self.config = self.make_config("GRCh37")

    def make_archive(self, assembly, entries=None):
        archive = self.sources / (assembly + "_v1.7.tar.gz")
        if entries is None:
            entries = [(assembly + "_v1.7/", None),
                       (assembly + "_v1.7/vep/", None),
                       (assembly + "_v1.7/vep/cache.txt", b"approved annotation\n"),
                       (assembly + "_v1.7/track.txt", b"unchanged track\n")]
        with tarfile.open(archive, "w:gz") as stream:
            for name, content in entries:
                info = tarfile.TarInfo(name)
                if content is None:
                    info.type = tarfile.DIRTYPE
                    stream.addfile(info)
                else:
                    info.size = len(content)
                    stream.addfile(info, io.BytesIO(content))
        return archive

    def make_config(self, assembly):
        archive = self.make_archive(assembly)
        prefix = "cadd_" + assembly + "_"
        config = {"assembly": assembly, "cadd_version": "v1.7",
                  "cache_root": str(self.directory / "cache"),
                  "cadd_parent_dir": str(self.base.parent),
                  prefix + "anno_url": archive.as_uri(),
                  prefix + "anno_md5": md5(archive.read_bytes())}
        for kind in ("snv", "indel"):
            data = self.sources / (assembly + "_" + kind + ".tsv.gz")
            data.write_bytes((kind + " approved data\n").encode())
            index = data.with_name(data.name + ".tbi")
            index.write_bytes((kind + " approved index\n").encode())
            config[prefix + kind + "_anno_url"] = data.as_uri()
            config[prefix + kind + "_anno_md5"] = md5(data.read_bytes())
            config[prefix + kind + "_anno_tbi_md5"] = md5(index.read_bytes())
        return config

    def install(self, config=None, assembly="GRCh37"):
        CADD.install(config or self.config, assembly, self.base)

    def annotation(self, assembly="GRCh37"):
        return self.base / "data/annotations" / (assembly + "_v1.7")

    def prescore(self, kind="snv", assembly="GRCh37"):
        return self.base / "data/prescored" / (assembly + "_v1.7") / "incl_anno" / (assembly + "_" + kind + ".tsv.gz")

    def test_complete_install_checks_data_and_each_index_and_retains_archive(self):
        self.install()
        self.assertTrue((self.annotation() / CADD.RECEIPT).is_file())
        self.assertTrue(self.annotation().with_name("GRCh37_v1.7.tar.gz").is_file())
        for kind in ("snv", "indel"):
            data = self.prescore(kind)
            for path in (data, data.with_name(data.name + ".tbi")):
                self.assertEqual(path.read_bytes(), (self.sources / path.name).read_bytes())
        self.assertEqual(list((self.base / "data").rglob("*.partial")), [])

    def test_index_missing_with_valid_data_recovers_independently(self):
        self.install()
        data = self.prescore()
        before = CADD.identity(data)
        index = data.with_name(data.name + ".tbi")
        index.unlink()
        with mock.patch.object(CADD, "urlopen", wraps=CADD.urlopen) as fetch:
            self.install()
        self.assertEqual(CADD.identity(data), before)
        self.assertEqual(index.read_bytes(), (self.sources / index.name).read_bytes())
        self.assertEqual([call.args[0].full_url for call in fetch.call_args_list],
                         [self.config["cadd_GRCh37_snv_anno_url"] + ".tbi"])

    def test_legacy_cluster_archive_reused_without_second_payload(self):
        legacy = Path(self.config["cache_root"]) / ".downloads/GRCh37_v1.7.tar.gz"
        legacy.parent.mkdir(parents=True)
        original = (self.sources / legacy.name).read_bytes()
        legacy.write_bytes(original)
        with mock.patch.object(CADD, "urlopen", wraps=CADD.urlopen) as fetch:
            self.install()
        self.assertFalse(self.annotation().with_name(legacy.name).exists())
        self.assertEqual(legacy.read_bytes(), original)
        receipt = json.loads((self.annotation() / CADD.RECEIPT).read_text())
        self.assertEqual(receipt["archive_identity"], CADD.receipt_identity(legacy))
        self.assertNotIn(self.config["cadd_GRCh37_anno_url"],
                         [call.args[0].full_url for call in fetch.call_args_list])

    def test_bad_canonical_archive_cannot_be_bypassed_by_legacy_copy(self):
        legacy = Path(self.config["cache_root"]) / ".downloads/GRCh37_v1.7.tar.gz"
        legacy.parent.mkdir(parents=True)
        legacy.write_bytes((self.sources / legacy.name).read_bytes())
        canonical = self.annotation().with_name(legacy.name)
        canonical.parent.mkdir(parents=True)
        canonical.write_bytes(b"unmatched retained archive")
        with mock.patch.object(CADD, "urlopen", side_effect=AssertionError("no replacement")), \
                self.assertRaisesRegex(ValueError, "checksum mismatch"):
            self.install()
        self.assertEqual(canonical.read_bytes(), b"unmatched retained archive")
        self.assertFalse(self.annotation().exists())

    def test_existing_bad_data_is_retained_and_install_fails(self):
        self.install()
        data = self.prescore()
        data.write_bytes(b"partial old final")
        with mock.patch.object(CADD, "urlopen", side_effect=AssertionError("must not overwrite")):
            with self.assertRaisesRegex(ValueError, "checksum mismatch.*existing bytes retained"):
                self.install()
        self.assertEqual(data.read_bytes(), b"partial old final")

    def test_existing_bad_index_is_retained_and_install_fails(self):
        self.install()
        index = self.prescore().with_name(self.prescore().name + ".tbi")
        index.write_bytes(b"wrong index")
        with self.assertRaisesRegex(ValueError, "checksum mismatch"):
            self.install()
        self.assertEqual(index.read_bytes(), b"wrong index")

    def test_missing_index_checksum_fails_before_any_cache_download(self):
        del self.config["cadd_GRCh37_snv_anno_tbi_md5"]
        with mock.patch.object(CADD, "urlopen", side_effect=AssertionError("must not download")):
            with self.assertRaisesRegex(ValueError, "cadd_GRCh37_snv_anno_tbi_md5"):
                self.install()
        self.assertFalse((self.base / "data").exists())

    def test_missing_annotation_url_is_explicit_and_never_retargets(self):
        del self.config["cadd_GRCh37_anno_url"]
        with self.assertRaisesRegex(ValueError, "Set cadd_GRCh37_anno_url"):
            self.install()

    def test_malformed_checksum_cannot_disable_validation(self):
        for value in (None, "", "null", ".*", "0" * 31):
            with self.subTest(value=value):
                config = dict(self.config, cadd_GRCh37_snv_anno_md5=value)
                with self.assertRaisesRegex(ValueError, "32-digit MD5"):
                    self.install(config)

    def test_index_checksum_is_not_the_data_checksum(self):
        self.config["cadd_GRCh37_snv_anno_tbi_md5"] = self.config["cadd_GRCh37_snv_anno_md5"]
        with self.assertRaisesRegex(ValueError, "checksum mismatch"):
            self.install()
        index = self.prescore().with_name(self.prescore().name + ".tbi")
        self.assertFalse(index.exists())
        self.assertTrue(index.with_name(index.name + ".partial").exists())
        self.assertFalse(self.prescore("indel").exists())

    def test_grch38_requires_its_own_fields(self):
        with self.assertRaisesRegex(ValueError, "cadd_GRCh38_anno_url"):
            self.install(assembly="hg38")
        config = self.make_config("GRCh38")
        self.install(config, assembly="hg38")
        self.assertTrue(self.prescore(assembly="GRCh38").is_file())
        self.assertFalse(self.prescore().exists())

    def test_unsupported_assembly_refused(self):
        with self.assertRaisesRegex(ValueError, "Unsupported CADD assembly"):
            self.install(assembly="hg18")

    def test_custom_urls_and_hashes_are_used_unchanged(self):
        resources = CADD.configured_resources(self.config, "hg19")
        self.assertEqual(resources["annotation"][0], self.config["cadd_GRCh37_anno_url"])
        self.assertEqual(resources["prescores"][0][0], self.config["cadd_GRCh37_snv_anno_url"])

    def test_query_url_requires_explicit_index_url(self):
        self.config["cadd_GRCh37_snv_anno_url"] += "?signed=fixture"
        with self.assertRaisesRegex(ValueError, "Set cadd_GRCh37_snv_anno_tbi_url explicitly"):
            self.install()
        self.config["cadd_GRCh37_snv_anno_tbi_url"] = (self.sources / "GRCh37_snv.tsv.gz.tbi").as_uri()
        self.assertEqual(CADD.configured_resources(self.config, "GRCh37")["prescores"][1][0],
                         self.config["cadd_GRCh37_snv_anno_tbi_url"])

    def test_completed_partial_is_published_without_an_http_request(self):
        destination = self.directory / "finished"
        partial = destination.with_name("finished.partial")
        partial.write_bytes(b"completed bytes")
        with mock.patch.object(CADD, "urlopen", side_effect=AssertionError("no request needed")):
            CADD.download("https://fixture.invalid/finished", destination, md5(partial.read_bytes()))
        self.assertEqual(destination.read_bytes(), b"completed bytes")
        self.assertFalse(partial.exists())

    def test_failed_download_never_publishes_final_name(self):
        destination = self.directory / "absent"
        with self.assertRaises(OSError):
            CADD.download((self.sources / "absent").as_uri(), destination, "0" * 32)
        self.assertFalse(destination.exists())

    def test_annotation_directory_without_resources_is_not_blessed(self):
        self.annotation().mkdir(parents=True)
        with mock.patch.object(CADD, "urlopen", side_effect=AssertionError("refuse early")):
            with self.assertRaisesRegex(ValueError, "Incomplete existing CADD annotations"):
                self.install()
        self.assertFalse((self.annotation() / CADD.RECEIPT).exists())

    def make_legacy_annotations(self):
        target = self.annotation()
        (target / "vep").mkdir(parents=True)
        (target / "vep/cache.txt").write_bytes(b"approved annotation\n")
        (target / "track.txt").write_bytes(b"unchanged track\n")
        return target

    def test_legacy_annotation_directory_is_fully_compared_then_reused(self):
        target = self.make_legacy_annotations()
        before = CADD.identity(target / "track.txt")
        self.install()
        self.assertEqual(CADD.identity(target / "track.txt"), before)
        self.assertTrue((target / CADD.RECEIPT).exists())

    def test_partial_legacy_tree_with_vep_is_refused_without_overwrite(self):
        target = self.make_legacy_annotations()
        (target / "track.txt").unlink()
        with self.assertRaisesRegex(ValueError, "existing CADD annotation"):
            self.install()
        self.assertFalse((target / "track.txt").exists())
        self.assertFalse((target / CADD.RECEIPT).exists())

    def test_changed_legacy_annotation_bytes_fail_full_comparison(self):
        target = self.make_legacy_annotations()
        (target / "track.txt").write_bytes(b"changed!! track\n")  # Same byte count.
        with self.assertRaisesRegex(ValueError, "differs from the verified archive"):
            self.install()
        self.assertEqual((target / "track.txt").read_bytes(), b"changed!! track\n")

    def test_unchanged_receipt_skips_rehashing_enormous_annotation_members(self):
        self.install()
        with mock.patch.object(CADD, "checksum", wraps=CADD.checksum) as hash_file:
            self.install()
        hashed = [Path(call.args[0]) for call in hash_file.call_args_list]
        self.assertFalse(any(path.is_relative_to(self.annotation()) for path in hashed))
        self.assertNotIn(self.annotation().with_name("GRCh37_v1.7.tar.gz"), hashed)
        self.assertEqual(len(hashed), 4)  # Existing prescore data and indexes still get MD5 gates.

    def node_metadata(self, device):
        """Simulate another client's device number without changing the object."""
        raw_identity, fstat = CADD.identity, CADD.os.fstat
        def raw(path):
            value = raw_identity(path)
            value[0] = device
            return value
        def descriptor(fd):
            value = fstat(fd)
            return SimpleNamespace(st_dev=device, **{name:getattr(value, name) for name in
                ("st_ino", "st_size", "st_mtime_ns", "st_ctime_ns", "st_mode")})
        stack = contextlib.ExitStack()
        stack.enter_context(mock.patch.object(CADD, "identity", side_effect=raw))
        stack.enter_context(mock.patch.object(CADD.os, "fstat", side_effect=descriptor))
        return stack

    def test_portable_reuse_preserves_receipt_across_device_and_node_changes(self):
        with mock.patch.object(CADD, "_filesystem_id", return_value=1234):
            self.install()
            receipt_path = self.annotation() / CADD.RECEIPT
            before = receipt_path.read_bytes(), CADD.identity(receipt_path)
            with self.node_metadata(987654), \
                    mock.patch.object(CADD.socket, "gethostname", return_value="other-node"), \
                    mock.patch.object(CADD, "checksum", wraps=CADD.checksum) as hashed, \
                    mock.patch.object(CADD.tempfile, "NamedTemporaryFile", side_effect=AssertionError("no republication")):
                self.install()
            self.assertEqual((receipt_path.read_bytes(), CADD.identity(receipt_path)), before)
            self.assertEqual(len(hashed.call_args_list), 4)
            self.assertTrue(all(self.prescore().parent == Path(call.args[0]).parent for call in hashed.call_args_list))

    def test_zero_or_unavailable_fsid_receipts_are_node_bound(self):
        path = self.sources / "GRCh37_v1.7.tar.gz"
        for metadata in (SimpleNamespace(f_fsid=0), SimpleNamespace()):
            with self.subTest(metadata=metadata), mock.patch.object(CADD.os, "fstatvfs", return_value=metadata):
                with mock.patch.object(CADD.socket, "gethostname", return_value="node-a"):
                    stored = CADD.receipt_identity(path)
                    self.assertEqual(stored, {"kind":"node-raw6-v1", "node":"node-a", "stat":CADD.identity(path)})
                    self.assertTrue(CADD.identity_matches(path, stored))
                with mock.patch.object(CADD.socket, "gethostname", return_value="node-b"):
                    self.assertFalse(CADD.identity_matches(path, stored))

    def test_unsupported_fsid_falls_back_but_io_errors_propagate(self):
        path = self.sources / "GRCh37_v1.7.tar.gz"
        with mock.patch.object(CADD.os, "fstatvfs", side_effect=OSError(errno.EOPNOTSUPP, "unsupported")):
            self.assertEqual(CADD.receipt_identity(path)["kind"], "node-raw6-v1")
        with mock.patch.object(CADD.os, "fstatvfs", side_effect=OSError(errno.EIO, "unavailable storage")):
            with self.assertRaises(OSError):
                CADD.receipt_identity(path)

    def test_symlink_identity_uses_its_own_filesystem_even_with_missing_target(self):
        link = self.directory / "link"
        link.symlink_to("missing-target")
        def filesystem(fd):
            self.assertTrue(CADD.stat.S_ISLNK(os.fstat(fd).st_mode))
            return 1234
        with mock.patch.object(CADD, "_filesystem_id", side_effect=filesystem) as captured:
            stored = CADD.receipt_identity(link)
        if hasattr(os, "O_PATH"):
            captured.assert_called_once()
            self.assertEqual(stored, {"kind":"statvfs-fsid-v1", "fsid":"1234", "stat":CADD.identity(link)[1:]})
        else:
            captured.assert_not_called()
            self.assertEqual(stored["kind"], "node-raw6-v1")
        self.assertTrue(link.is_symlink())
        self.assertFalse(link.exists())

    def test_unknown_malformed_and_changed_identities_cannot_skip_archive_md5(self):
        path = self.sources / "GRCh37_v1.7.tar.gz"
        with mock.patch.object(CADD, "_filesystem_id", return_value=1234):
            known = CADD.receipt_identity(path)
            malformed = [None, {}, {"kind":"future-v2"},
                         dict(known, fsid="0"), dict(known, fsid="01234"), dict(known, fsid=1234),
                         dict(known, stat=[True]+known["stat"][1:]), dict(known, extra="unknown"),
                         {"kind":"node-raw6-v1", "node":True, "stat":CADD.identity(path)},
                         [True]+CADD.identity(path)[1:]]
            for index in range(5):
                changed = list(known["stat"]); changed[index] += 1
                malformed.append(dict(known, stat=changed))
            malformed.append(dict(known, fsid="5678"))
            for stored in malformed:
                with self.subTest(stored=stored), mock.patch.object(CADD, "checksum", wraps=CADD.checksum) as hashed:
                    with self.assertRaisesRegex(ValueError, "checksum mismatch"):
                        CADD.download("https://fixture.invalid/archive", path, "0"*32, stored)
                    hashed.assert_called_once_with(path)

    def test_legacy_raw6_migration_requires_the_complete_tuple(self):
        self.install()
        receipt_path = self.annotation() / CADD.RECEIPT
        archive = self.annotation().with_name("GRCh37_v1.7.tar.gz")
        asset = self.annotation() / "track.txt"
        def legacy(device_mismatch):
            record = json.loads(receipt_path.read_text())
            record["archive_identity"] = CADD.identity(archive)
            for name, entry in record["members"].items():
                if entry["kind"] == "file":
                    entry["identity"] = CADD.identity(self.annotation()/name)
            if device_mismatch:
                record["archive_identity"][0] += 1
                record["members"]["track.txt"]["identity"][0] += 1
            receipt_path.write_text(json.dumps(record))
        legacy(False)
        with mock.patch.object(CADD, "checksum", wraps=CADD.checksum) as hashed:
            self.install()
        self.assertEqual(len(hashed.call_args_list), 4)
        self.assertIsInstance(json.loads(receipt_path.read_text())["archive_identity"], dict)
        legacy(True)
        with mock.patch.object(CADD, "checksum", wraps=CADD.checksum) as hashed:
            self.install()
        self.assertIn(mock.call(archive), hashed.call_args_list)
        self.assertIn(mock.call(asset, "sha256"), hashed.call_args_list)
        self.assertEqual(len(hashed.call_args_list), 6)

    def test_fsid_change_reverifies_archive_and_updates_receipt(self):
        with mock.patch.object(CADD, "_filesystem_id", return_value=1234):
            self.install()
        with mock.patch.object(CADD, "_filesystem_id", return_value=5678), \
                mock.patch.object(CADD, "checksum", wraps=CADD.checksum) as hashed:
            self.install()
        self.assertIn(mock.call(self.annotation().with_name("GRCh37_v1.7.tar.gz")), hashed.call_args_list)
        receipt = json.loads((self.annotation()/CADD.RECEIPT).read_text())
        self.assertEqual(receipt["archive_identity"]["fsid"], "5678")
        self.assertEqual(receipt["members"]["track.txt"]["identity"]["fsid"], "5678")

    def test_identity_race_fails_without_blessing_a_local_fallback(self):
        path = self.sources / "GRCh37_v1.7.tar.gz"
        def race(_fd):
            value = path.stat()
            os.utime(path, ns=(value.st_atime_ns, value.st_mtime_ns+1_000_000_000))
            return None
        with mock.patch.object(CADD, "_filesystem_id", side_effect=race), \
                self.assertRaisesRegex(ValueError, "changed during identity verification"):
            CADD.receipt_identity(path)

    def test_path_replacement_during_identity_capture_is_refused(self):
        path = self.sources / "GRCh37_v1.7.tar.gz"
        replacement = path.with_name("replacement")
        replacement.write_bytes(path.read_bytes())
        def race(_fd):
            replacement.replace(path)
            return 1234
        with mock.patch.object(CADD, "_filesystem_id", side_effect=race), \
                self.assertRaisesRegex(ValueError, "changed during identity verification"):
            CADD.receipt_identity(path)

    def test_member_metadata_race_preserves_existing_ready_receipt(self):
        self.install()
        receipt_path = self.annotation()/CADD.RECEIPT
        before = receipt_path.read_bytes(), CADD.identity(receipt_path)
        asset = self.annotation()/"track.txt"
        inode = asset.stat().st_ino
        original = CADD._filesystem_id
        def race(fd):
            value = original(fd)
            if os.fstat(fd).st_ino == inode:
                current = asset.stat()
                os.utime(asset, ns=(current.st_atime_ns, current.st_mtime_ns+1_000_000_000))
            return value
        with mock.patch.object(CADD, "_filesystem_id", side_effect=race), \
                self.assertRaisesRegex(ValueError, "changed during identity verification"):
            self.install()
        self.assertEqual((receipt_path.read_bytes(), CADD.identity(receipt_path)), before)

    def test_verified_reuse_preserves_receipt_and_annotation_identities_for_both_assemblies(self):
        for assembly in ("GRCh37", "GRCh38"):
            with self.subTest(assembly=assembly):
                config = self.config if assembly == "GRCh37" else self.make_config(assembly)
                self.install(config, assembly)
                target = self.annotation(assembly)
                paths = [target, *target.rglob("*")]
                before = {path: CADD.identity(path) for path in paths}
                receipt_bytes = (target / CADD.RECEIPT).read_bytes()
                with mock.patch.object(CADD, "verify_receipt", wraps=CADD.verify_receipt) as verify, \
                        mock.patch.object(CADD.tempfile, "NamedTemporaryFile", side_effect=AssertionError("unchanged receipt must not be republished")), \
                        mock.patch.object(CADD, "urlopen", side_effect=AssertionError("no download")):
                    self.install(config, assembly)
                verify.assert_called_once()
                self.assertEqual(before, {path: CADD.identity(path) for path in paths})
                self.assertEqual((target / CADD.RECEIPT).read_bytes(), receipt_bytes)
                self.assertEqual(set(target.rglob("*")), set(paths[1:]))

    def test_verified_archive_identity_change_updates_receipt_then_reuses_it(self):
        self.install()
        archive = self.annotation().with_name("GRCh37_v1.7.tar.gz")
        receipt_path = self.annotation() / CADD.RECEIPT
        before = receipt_path.read_bytes()
        stat = archive.stat()
        os.utime(archive, ns=(stat.st_atime_ns, stat.st_mtime_ns + 1_000_000_000))
        with mock.patch.object(CADD, "checksum", wraps=CADD.checksum) as hash_file:
            self.install()
        self.assertIn(archive, [Path(call.args[0]) for call in hash_file.call_args_list])
        self.assertNotEqual(receipt_path.read_bytes(), before)
        self.assertEqual(json.loads(receipt_path.read_bytes())["archive_identity"], CADD.receipt_identity(archive))
        published = CADD.identity(receipt_path), receipt_path.read_bytes()
        self.install()
        self.assertEqual((CADD.identity(receipt_path), receipt_path.read_bytes()), published)

    def test_changed_identity_rechecks_content_and_refuses_changed_annotations(self):
        self.install()
        target = self.annotation() / "track.txt"
        target.write_bytes(b"changed!! track\n")
        with self.assertRaisesRegex(ValueError, "Existing CADD annotation changed"):
            self.install()
        self.assertEqual(target.read_bytes(), b"changed!! track\n")

    def test_identical_bytes_with_changed_identity_are_verified_and_receipt_updated(self):
        self.install()
        target = self.annotation() / "track.txt"
        receipt_path = self.annotation() / CADD.RECEIPT
        before = receipt_path.read_bytes()
        stat = target.stat()
        os.utime(target, ns=(stat.st_atime_ns, stat.st_mtime_ns + 1_000_000_000))
        with mock.patch.object(CADD, "checksum", wraps=CADD.checksum) as hash_file:
            self.install()
        self.assertIn(target, [Path(call.args[0]) for call in hash_file.call_args_list])
        receipt = json.loads((self.annotation() / CADD.RECEIPT).read_text())
        self.assertEqual(receipt["members"]["track.txt"]["identity"], CADD.receipt_identity(target))
        self.assertNotEqual(receipt_path.read_bytes(), before)
        published = CADD.identity(receipt_path), receipt_path.read_bytes()
        self.install()
        self.assertEqual((CADD.identity(receipt_path), receipt_path.read_bytes()), published)

    def test_missing_receipt_member_is_refused(self):
        self.install()
        target = self.annotation() / "track.txt"
        target.unlink()
        with self.assertRaisesRegex(ValueError, "Incomplete existing CADD annotations"):
            self.install()
        self.assertFalse(target.exists())

    def test_bad_archive_does_not_publish_extracted_directory(self):
        archive = self.sources / "bad.tar.gz"
        archive.write_bytes(b"not a tar")
        self.config["cadd_GRCh37_anno_url"] = archive.as_uri()
        self.config["cadd_GRCh37_anno_md5"] = md5(archive.read_bytes())
        with self.assertRaises(tarfile.TarError):
            self.install()
        self.assertFalse(self.annotation().exists())
        self.assertEqual(list(self.annotation().parent.glob(".cadd-extract-*")), [])

    def test_archive_metadata_change_after_real_extraction_removes_staging(self):
        extract = CADD.verify_or_extract
        def extract_then_change_archive(archive, target, extracting):
            members = extract(archive, target, extracting)
            self.assertTrue(extracting)
            self.assertTrue((target/"vep/cache.txt").is_file())
            value = archive.stat()
            os.utime(archive, ns=(value.st_atime_ns, value.st_mtime_ns+1_000_000_000))
            return members
        with mock.patch.object(CADD, "verify_or_extract", side_effect=extract_then_change_archive), \
                self.assertRaisesRegex(ValueError, "CADD cache changed during verification"):
            self.install()
        self.assertFalse(self.annotation().exists())
        self.assertEqual(list(self.annotation().parent.glob(".cadd-extract-*")), [])

    def test_archive_path_traversal_refused_in_staging(self):
        archive = self.make_archive("GRCh37", [("GRCh37_v1.7/../escaped.txt", b"escape")])
        self.config["cadd_GRCh37_anno_md5"] = md5(archive.read_bytes())
        with self.assertRaisesRegex(ValueError, "unexpected member path"):
            self.install()
        self.assertFalse(self.annotation().exists())
        self.assertFalse((self.annotation().parent / "escaped.txt").exists())

    def link_archive(self, link_target):
        archive = self.sources / "GRCh37_v1.7.tar.gz"
        with tarfile.open(archive, "w:gz") as stream:
            content = b"approved linked annotation\n"
            regular = tarfile.TarInfo("GRCh37_v1.7/vep/track.txt")
            regular.size = len(content)
            stream.addfile(regular, io.BytesIO(content))
            symlink = tarfile.TarInfo("GRCh37_v1.7/vep/symlink.txt")
            symlink.type, symlink.linkname = tarfile.SYMTYPE, link_target
            stream.addfile(symlink)
            hardlink = tarfile.TarInfo("GRCh37_v1.7/hardlink.txt")
            hardlink.type, hardlink.linkname = tarfile.LNKTYPE, "GRCh37_v1.7/vep/track.txt"
            stream.addfile(hardlink)
        self.config["cadd_GRCh37_anno_md5"] = md5(archive.read_bytes())

    def test_safe_archive_links_survive_publication_and_receipt_reuse(self):
        self.link_archive("track.txt")
        self.install()
        target = self.annotation()
        self.assertEqual((target / "vep/symlink.txt").read_bytes(), (target / "hardlink.txt").read_bytes())
        with mock.patch.object(CADD, "checksum", wraps=CADD.checksum) as hash_file:
            self.install()
        self.assertFalse(any(Path(call.args[0]).is_relative_to(target) for call in hash_file.call_args_list))

    def test_archive_symlink_cannot_escape_annotation_directory(self):
        self.link_archive("../../escaped.txt")
        with self.assertRaisesRegex(ValueError, "symlink escapes"):
            self.install()
        self.assertFalse(self.annotation().exists())
        self.assertFalse((self.annotation().parent / "escaped.txt").exists())

    def shell_install(self, assembly="hg19", legacy=False, configure=False, legacy_parent=None):
        config = dict(self.config)
        if legacy:
            config.pop("cache_root")
        path = self.directory / "config.json"
        path.write_text(json.dumps(config))
        script = r'''
            source "$1" >/dev/null || exit
            python() { "$PRIVA_TEST_PYTHON" "$@"; }
            if [[ "$PRIVA_TEST_CONFIGURE" != yes ]]; then
                priva_configure_cache_paths() { return 0; }
            fi
            update_or_append_yaml() {
                python -c 'import json,sys; from pathlib import Path; path=Path(sys.argv[1]); data=json.loads(path.read_text()); data[sys.argv[2]]=sys.argv[3]; path.write_text(json.dumps(data))' "$@"
            }
            yq() { python -c 'import json,sys; print(json.dumps(json.load(open(sys.argv[1]))))' "$2"; }
            read_yaml() { python -c 'import json,sys; print(json.load(open(sys.argv[1])).get(sys.argv[2], ""))' "$1" "$2"; }
            if [[ "$PRIVA_TEST_LEGACY" == yes ]]; then
                read() {
                    if [[ "$4" == reply ]]; then
                        printf -v "$4" '%s' "$PRIVA_TEST_LEGACY_REPLY"
                    else
                        printf -v "$4" '%s' "$PRIVA_TEST_LEGACY_PARENT"
                    fi
                }
            else
                read() { printf 'UNEXPECTED_PROMPT\n' >&2; return 90; }
            fi
            CADD_install "$2" "$3"
        '''
        return subprocess.run(["bash", "--noprofile", "--norc", "-c", script, "bash",
                               str(ROOT / "scripts/install_utils.sh"), str(path), assembly],
                              env=dict(os.environ, PRIVA_TEST_PYTHON=sys.executable,
                                       PRIVA_TEST_LEGACY="yes" if legacy else "no",
                                       PRIVA_TEST_CONFIGURE="yes" if configure else "no",
                                       PRIVA_TEST_LEGACY_REPLY="no" if legacy_parent else "yes",
                                       PRIVA_TEST_LEGACY_PARENT=str(legacy_parent or ""),
                                       CONDA_PREFIX=str(self.base)),
                              text=True, capture_output=True, timeout=30)

    def test_shell_configured_parent_routes_to_helper_without_prompt(self):
        result = self.shell_install()
        self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
        self.assertNotIn("UNEXPECTED_PROMPT", result.stdout + result.stderr)
        self.assertTrue(self.prescore().is_file())

    def test_shell_selected_nested_source_precedes_matching_sibling_without_prompt(self):
        selected = self.base.parent / "nested" / "custom checkout"
        selected.mkdir(parents=True); (selected / "CADD.sh").write_bytes((self.base / "CADD.sh").read_bytes())
        self.config["cadd_base_dir"] = str(selected)
        before = CADD.identity(selected / "CADD.sh")
        result = self.shell_install()
        self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
        self.assertEqual(json.loads((self.directory / "config.json").read_text())["cadd_base_dir"], str(selected))
        self.assertTrue((selected / "data/annotations/GRCh37_v1.7/track.txt").is_file())
        self.assertFalse((self.base / "data").exists())
        self.assertEqual(CADD.identity(selected / "CADD.sh"), before)

    def test_shell_planner_adopted_supplied_source_stays_selected(self):
        source = self.sources / "selected-source"; source.mkdir()
        (source / "CADD.sh").write_bytes((self.base / "CADD.sh").read_bytes())
        (source / "original-marker").write_bytes(b"supplied source identity\n")
        before = {path.name: (CADD.identity(path), path.read_bytes()) for path in source.iterdir()}
        self.config.pop("cadd_parent_dir")
        self.config.update(tmp_dir=str(self.directory / "tmp"), vep_installed_version="116")
        self.config["cadd_base_dir"] = str(source)
        parent = Path(self.config["cache_root"]) / "CADD"
        canonical = parent / "CADD-scripts-1.7.2"; canonical.mkdir(parents=True)
        (canonical / "CADD.sh").write_text("canonical decoy must not be selected\n")
        result = self.shell_install(configure=True)
        self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
        selected = parent / "supplied/cadd_base_dir/selected-source"
        self.assertEqual(json.loads((self.directory / "config.json").read_text())["cadd_base_dir"], str(selected))
        self.assertTrue((selected / "data/annotations/GRCh37_v1.7/track.txt").is_file())
        self.assertEqual((selected / "original-marker").read_bytes(), before["original-marker"][1])
        self.assertFalse((canonical / "data").exists())
        self.assertEqual(before, {path.name: (CADD.identity(path), path.read_bytes()) for path in source.iterdir()})

    def test_shell_invalid_selected_source_is_retained_without_using_matching_sibling(self):
        empty = self.directory / "empty-source"; empty.mkdir()
        regular = self.directory / "source-file"; regular.write_bytes(b"retain this selected file\n")
        missing = self.base.parent / "custom-missing-source"
        broken = self.directory / "broken-source"; broken.symlink_to(missing)
        for selected in (missing, empty, regular, broken):
            with self.subTest(selected=selected.name):
                self.config["cadd_base_dir"] = str(selected)
                result = self.shell_install()
                self.assertNotEqual(result.returncode, 0)
                self.assertIn("selected path retained: " + str(selected), result.stdout + result.stderr)
                self.assertEqual(json.loads((self.directory / "config.json").read_text()), self.config)
                self.assertFalse((self.base / "data").exists())
        self.assertEqual(regular.read_bytes(), b"retain this selected file\n")
        self.assertTrue(broken.is_symlink()); self.assertEqual(list(empty.iterdir()), [])

    def test_shell_helper_failure_propagates_even_without_errexit(self):
        self.config["cadd_GRCh37_anno_md5"] = "0" * 32
        result = self.shell_install()
        self.assertNotEqual(result.returncode, 0)
        self.assertIn("CADD cache installation failed", result.stderr)
        self.assertFalse(self.annotation().exists())
        self.assertFalse(self.prescore().exists())

    def test_shell_legacy_environment_selection_is_preserved(self):
        result = self.shell_install(legacy=True)
        self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
        self.assertTrue(self.prescore().exists())

    def test_shell_legacy_alternate_parent_overrides_only_environment_discovery(self):
        parent = self.directory / "legacy alternate parent"; selected = parent / self.base.name
        selected.mkdir(parents=True); (selected / "CADD.sh").write_bytes((self.base / "CADD.sh").read_bytes())
        result = self.shell_install(legacy=True, legacy_parent=parent)
        self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
        self.assertEqual(json.loads((self.directory / "config.json").read_text())["cadd_base_dir"], str(selected))
        self.assertTrue((selected / "data/annotations/GRCh37_v1.7/track.txt").is_file())
        self.assertFalse((self.base / "data").exists())

    def test_shell_fresh_source_uses_reviewed_git_instead_of_unverified_zip(self):
        (self.base / "CADD.sh").unlink()
        self.base.rmdir()
        self.config["cadd_base_dir"] = str(self.base)
        self.config["cadd_zip_download_url"] = "https://github.com/kircherlab/CADD-scripts/archive/refs/tags/v1.7.2.zip"
        binary = self.directory / "fake-git-bin"
        binary.mkdir()
        git = binary / "git"
        git.write_text("#!" + sys.executable + "\n" + '''import json, os, pathlib, sys
args = sys.argv[1:]
with open(os.environ['CADD_TEST_GIT_LOG'], 'a') as stream:
    stream.write(json.dumps(args) + '\\n')
if args[0] == 'clone':
    destination = pathlib.Path(args[-1])
    destination.mkdir()
    (destination / '.git').mkdir()
    (destination / 'CADD.sh').write_text('#!/bin/bash\\nexit 99\\n')
elif args[0] == '-C' and args[2:] == ['rev-parse', 'HEAD']:
    print('1fead0c1e02d3f0a8599bdb04e6b3a2dd70b4866')
elif args[0] != '-C' or args[2] not in ('checkout', 'diff', 'ls-files'):
    raise SystemExit('Unexpected fake Git invocation: ' + repr(args))
''')
        git.chmod(0o755)
        log = self.directory / "git-commands.jsonl"
        with mock.patch.dict(os.environ, {"PATH": str(binary) + os.pathsep + os.environ["PATH"],
                                          "CADD_TEST_GIT_LOG": str(log)}):
            result = self.shell_install()
        self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
        commands = [json.loads(line) for line in log.read_text().splitlines()]
        self.assertEqual(commands[0], ['clone', '--no-checkout', 'https://github.com/kircherlab/CADD-scripts.git', str(self.base)])
        self.assertIn(['-C', str(self.base), 'checkout', '--detach', '1fead0c1e02d3f0a8599bdb04e6b3a2dd70b4866'], commands)
        self.assertFalse((self.base.parent / 'CADD-scripts.zip').exists())
        self.assertTrue(self.prescore().exists())

    def test_shell_fresh_custom_archive_fails_without_replacing_source(self):
        (self.base / "CADD.sh").unlink()
        self.base.rmdir()
        self.config["cadd_zip_download_url"] = "https://example.invalid/custom-cadd.zip"
        result = self.shell_install()
        self.assertNotEqual(result.returncode, 0)
        self.assertIn("not the reviewed v1.7.2 source", result.stderr)
        self.assertFalse(self.base.exists())


class CaddHttpRecoveryTests(unittest.TestCase):
    def setUp(self):
        temporary = tempfile.TemporaryDirectory(prefix="priva_cadd_http_", dir=os.environ.get("PRIVA_TEST_TMPDIR", os.environ.get("TMPDIR")))
        self.addCleanup(temporary.cleanup)
        self.directory = Path(temporary.name)
        self.payload = b"tiny approved HTTP payload\n"
        self.requests = []
        self.mode = "range"
        owner = self

        class Handler(BaseHTTPRequestHandler):
            def do_GET(self):
                request_range = self.headers.get("Range")
                owner.requests.append(request_range)
                offset = int(request_range[6:-1]) if request_range else 0
                if offset >= len(owner.payload):
                    self.send_error(416)
                elif owner.mode == "short":
                    self.send_response(200)
                    self.send_header("Content-Length", str(len(owner.payload)))
                    self.end_headers()
                    self.wfile.write(owner.payload[:7])
                else:
                    resumed = request_range and owner.mode == "range"
                    content = owner.payload[offset:] if resumed else owner.payload
                    self.send_response(206 if resumed else 200)
                    if resumed:
                        self.send_header("Content-Range", f"bytes {offset}-{len(owner.payload) - 1}/{len(owner.payload)}")
                    self.send_header("Content-Length", str(len(content)))
                    self.end_headers()
                    self.wfile.write(content)

            def log_message(self, *args):
                pass

        self.server = ThreadingHTTPServer(("127.0.0.1", 0), Handler)
        self.thread = threading.Thread(target=self.server.serve_forever, daemon=True)
        self.thread.start()
        self.addCleanup(self.stop_server)
        self.url = f"http://127.0.0.1:{self.server.server_port}/payload"
        self.destination = self.directory / "payload"
        self.partial = self.directory / "payload.partial"

    def stop_server(self):
        self.server.shutdown()
        self.server.server_close()
        self.thread.join(timeout=5)

    def download(self, checksum=None):
        with contextlib.redirect_stdout(io.StringIO()):
            CADD.download(self.url, self.destination, checksum or md5(self.payload))

    def test_prefix_partial_resumes_then_publishes_verified_final(self):
        self.partial.write_bytes(self.payload[:7])
        self.download()
        self.assertEqual(self.requests, ["bytes=7-"])
        self.assertEqual(self.destination.read_bytes(), self.payload)
        self.assertFalse(self.partial.exists())

    def test_server_ignoring_range_restarts_partial_safely(self):
        self.mode = "ignore"
        self.partial.write_bytes(b"bad old prefix")
        self.download()
        self.assertEqual(self.destination.read_bytes(), self.payload)

    def test_checksum_failure_leaves_partial_and_no_final(self):
        with self.assertRaisesRegex(ValueError, "checksum mismatch"):
            self.download("0" * 32)
        self.assertFalse(self.destination.exists())
        self.assertEqual(self.partial.read_bytes(), self.payload)

    def test_short_response_is_not_published_and_rerun_recovers(self):
        self.mode = "short"
        with self.assertRaises((OSError, ValueError, CADD.http.client.HTTPException)):
            self.download()
        self.assertFalse(self.destination.exists())
        self.assertTrue(self.partial.exists())
        self.assertEqual(self.partial.read_bytes(), self.payload[:7])
        self.mode = "range"
        self.download()
        self.assertEqual(self.destination.read_bytes(), self.payload)

    def test_invalid_complete_partial_is_retained_with_actionable_resume_error(self):
        self.partial.write_bytes(b"x" * len(self.payload))
        with self.assertRaisesRegex(ValueError, "move the partial aside explicitly"):
            self.download()
        self.assertFalse(self.destination.exists())
        self.assertEqual(self.partial.read_bytes(), b"x" * len(self.payload))


if __name__ == "__main__":
    unittest.main()
