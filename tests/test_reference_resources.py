"""Pinned reference preparation with tiny real FASTAs and mocked HTTP bytes."""
import copy
import contextlib
import gzip
import hashlib
import importlib.util
import io
import json
import os
from pathlib import Path
import pickle
import shutil
import subprocess
import sys
import tempfile
import unittest
from unittest import mock

import fake_priva_binary


ROOT = Path(__file__).resolve().parents[1]
SPEC = importlib.util.spec_from_file_location("reference_resources", ROOT / "scripts/reference_resources.py")
REF = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(REF)
PRE_SPEC = importlib.util.spec_from_file_location("reference_prerequisites", ROOT / "scripts/installer_prerequisites.py")
PRE = importlib.util.module_from_spec(PRE_SPEC)
PRE_SPEC.loader.exec_module(PRE)


class Response(io.BytesIO):
    def __init__(self, body, status=200, headers=None):
        super().__init__(body)
        self.status = status
        self.headers = headers or {}

    def getcode(self):
        return self.status


class InterruptedResponse(Response):
    def read(self, amount=-1):
        if self.tell():
            raise OSError("fixture connection interrupted")
        return super().read(8)


class ReferenceResources(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory(prefix="reference_resources_",
                                               dir=os.environ.get("PRIVA_TEST_TMPDIR", os.environ.get("TMPDIR")))
        self.addCleanup(self.temp.cleanup)
        self.folder = Path(self.temp.name)
        self.data = self.folder / "selected_source/data"
        self.patch_lengths = mock.patch.dict(REF.CHR1_LENGTH, {"hg19": 12, "hg38": 14})
        self.patch_lengths.start()
        self.addCleanup(self.patch_lengths.stop)
        self.policy = mock.patch.dict(os.environ, {"PRIVA_INSTALL_INSPECT_ROOTS": str(self.folder)})
        self.policy.start()
        self.addCleanup(self.policy.stop)

    def reference(self, assembly="hg38", name=None, index=True):
        path = self.folder / (name or f"original.{assembly}.fasta")
        path.write_text(">chr1\n" + "A" * REF.CHR1_LENGTH[assembly] + "\n>chr2\nACGT\n")
        if index:
            if not shutil.which("samtools"):
                self.skipTest("samtools is required for tiny real FASTA fixtures")
            subprocess.run(["samtools", "faidx", str(path)], check=True, capture_output=True)
        return path

    def sources(self):
        bodies = {
            "hg19": b">chr1\nAAAAAAAAAAAA\n>chr2\nACGT\n",
            "hg38": b">chr1\nAAAAAAAAAAAAAA\n>chr2\nACGT\n",
            "hg19_to_hg38": b"chain 1 chr1 12 + 0 1 chr1 14 + 0 1 1\n1\n",
        }
        sources = copy.deepcopy(REF.SOURCES)
        archives = {}
        for key, body in bodies.items():
            compressed = gzip.compress(body, mtime=0)
            archives[sources[key]["url"]] = compressed
            sources[key]["md5"] = hashlib.md5(compressed).hexdigest()
            sources[key]["bytes"] = len(compressed)
        return sources, bodies, archives

    def installed(self):
        if not shutil.which("samtools"):
            self.skipTest("samtools is required for tiny real FASTA fixtures")
        sources, bodies, archives = self.sources()
        resources, updates = REF.resource_plan({"assembly": "hg38"}, self.data)
        patch = mock.patch.dict(REF.SOURCES, sources)
        patch.start()
        self.addCleanup(patch.stop)
        with mock.patch.object(REF, "urlopen", side_effect=lambda request, **_: Response(archives[request.full_url])) as network:
            REF.preflight(resources)
            for resource in resources:
                REF.install_resource(resource)
        self.assertEqual(network.call_count, 3)
        return resources, updates, bodies

    def test_blank_defaults_resolve_both_assemblies_without_writes(self):
        resources, updates = REF.resource_plan({"assembly": "GRCh38"}, self.data)
        REF.preflight(resources)
        self.assertEqual(len(resources), 3)
        self.assertEqual(updates["ref_genome"], updates["gofcards_hg38_fasta"])
        self.assertNotEqual(updates["gofcards_hg19_fasta"], updates["gofcards_hg38_fasta"])
        self.assertFalse(self.data.exists())

    def test_filename_never_selects_the_assembly(self):
        with self.assertRaisesRegex(ValueError, "filenames do not select"):
            REF.resource_plan({"ref_genome": "ucsc.hg19.fasta"}, self.data)

    def test_opt_in_cache_root_routes_fresh_defaults_for_all_assembly_aliases(self):
        cache = self.folder / "new_cache"
        for assembly, canonical in (("hg19", "hg19"), ("GRCh37", "hg19"),
                                    ("hg38", "hg38"), ("GRCh38", "hg38")):
            with self.subTest(assembly=assembly):
                config = {"assembly": assembly, "cache_root": str(cache)}
                original = copy.deepcopy(config)
                resources, updates = REF.resource_plan(config, self.data, use_cache_root=True)
                REF.preflight(resources)
                self.assertEqual(updates["ref_genome"], str(cache / "reference" / f"ucsc.{canonical}.fasta"))
                self.assertEqual(updates[f"gofcards_{canonical}_fasta"], updates["ref_genome"])
                self.assertEqual(updates["gofcards_hg19_to_hg38_chain"], str(cache / "liftOver/hg19ToHg38.over.chain"))
                self.assertTrue(all(resource["managed"] for resource in resources))
                self.assertEqual(config, original)
                self.assertFalse(cache.exists())
                self.assertFalse(self.data.exists())

    def test_direct_data_dir_semantics_are_unchanged_without_opt_in(self):
        cache = self.folder / "new_cache"
        config = {"assembly": "hg38", "cache_root": str(cache)}
        resources, updates = REF.resource_plan(config, self.data, "selected")
        self.assertEqual(updates["ref_genome"], str(self.data / "reference/ucsc.hg38.fasta"))
        config_path = self.folder / "config.yaml"
        config_path.write_text(f"assembly: hg38\ncache_root: {cache}\n")
        for opt_in, expected in ((False, self.data), (True, cache)):
            with self.subTest(opt_in=opt_in):
                argv = ["reference_resources.py", "preflight", str(config_path), "--data-dir", str(self.data),
                        "--scope", "selected"] + (["--use-cache-root"] if opt_in else [])
                with mock.patch.object(sys, "argv", argv), mock.patch.object(REF, "preflight") as preflight:
                    self.assertEqual(REF.main(), 0)
                self.assertEqual(preflight.call_args.args[0][0]["path"], expected / "reference/ucsc.hg38.fasta")
                self.assertFalse(cache.exists())
                self.assertFalse(self.data.exists())

    def test_only_absent_exact_legacy_defaults_redirect_to_cache_root(self):
        cache = self.folder / "new_cache"
        config = {"assembly": "hg19", "cache_root": str(cache),
                  "ref_genome": str(self.data / "reference/ucsc.hg19.fasta"),
                  "gofcards_hg38_fasta": str(self.data / "reference/ucsc.hg38.fasta"),
                  "gofcards_hg19_to_hg38_chain": str(self.data / "liftOver/hg19ToHg38.over.chain")}
        before = copy.deepcopy(config)
        resources, updates = REF.resource_plan(config, self.data, use_cache_root=True)
        REF.preflight(resources)
        self.assertTrue(all(cache in resource["path"].parents for resource in resources))
        self.assertTrue(all(resource["managed"] for resource in resources))
        self.assertEqual(config, before)
        self.assertEqual(updates["gofcards_hg19_fasta"], updates["ref_genome"])
        self.assertFalse(cache.exists())
        self.assertFalse(self.data.exists())

    def test_explicit_new_canonical_paths_are_managed_without_filename_guessing(self):
        cache = self.folder / "new_cache"
        config = {"assembly": "GRCh38", "cache_root": str(cache),
                  "ref_genome": str(cache / "reference/ucsc.hg38.fasta"),
                  "gofcards_hg19_to_hg38_chain": str(cache / "liftOver/hg19ToHg38.over.chain")}
        resources, updates = REF.resource_plan(config, self.data, use_cache_root=True)
        REF.preflight(resources)
        self.assertTrue(all(resource["managed"] for resource in resources))
        self.assertEqual(updates["ref_genome"], config["ref_genome"])
        self.assertFalse(cache.exists())

    def test_existing_original_selection_is_preserved_with_new_opposite_defaults(self):
        original = self.reference(name="selected_original.fa")
        before = (original.read_bytes(), Path(str(original) + ".fai").read_bytes())
        cache = self.folder / "new_cache"
        config = {"assembly": "hg38", "cache_root": str(cache), "ref_genome": str(original)}
        resources, updates = REF.resource_plan(config, self.data, use_cache_root=True)
        REF.preflight(resources)
        selected = next(resource for resource in resources if resource["path"] == original)
        self.assertFalse(selected["managed"])
        self.assertEqual(updates["ref_genome"], str(original))
        self.assertEqual(updates["gofcards_hg38_fasta"], str(original))
        self.assertEqual(updates["gofcards_hg19_fasta"], str(cache / "reference/ucsc.hg19.fasta"))
        self.assertEqual((original.read_bytes(), Path(str(original) + ".fai").read_bytes()), before)
        self.assertFalse(cache.exists())

    def test_existing_legacy_defaults_retain_managed_source_and_receipt_checks(self):
        previous, _, _ = self.installed()
        cache = self.folder / "new_cache"
        config = {"assembly": "hg38", "cache_root": str(cache),
                  "ref_genome": str(previous[0]["path"]),
                  "gofcards_hg19_fasta": str(previous[1]["path"]),
                  "gofcards_hg19_to_hg38_chain": str(previous[2]["path"])}
        resources, updates = REF.resource_plan(config, self.data, use_cache_root=True)
        self.assertEqual([resource["path"] for resource in resources], [resource["path"] for resource in previous])
        self.assertTrue(all(resource["managed"] for resource in resources))
        with mock.patch.object(REF, "verify_managed", wraps=REF.verify_managed) as verified:
            REF.preflight(resources)
        self.assertEqual(verified.call_count, 3)
        self.assertEqual(updates["ref_genome"], config["ref_genome"])
        self.assertFalse(cache.exists())

    def test_cache_root_never_downgrades_corrupt_old_managed_reference_to_custom(self):
        previous, _, _ = self.installed()
        old = previous[0]["path"]
        altered = old.read_bytes().replace(b"A", b"C", 1)
        old.write_bytes(altered)
        cache = self.folder / "new_cache"
        resources, _ = REF.resource_plan({"assembly": "hg38", "cache_root": str(cache), "ref_genome": str(old)},
                                        self.data, "selected", use_cache_root=True)
        with self.assertRaisesRegex(ValueError, "reference checksum mismatch"):
            REF.preflight(resources)
        self.assertEqual(old.read_bytes(), altered)
        self.assertFalse(cache.exists())

    def test_existing_unproven_old_managed_reference_is_not_relabelled_as_custom(self):
        old = self.data / "reference/ucsc.hg38.fasta"
        old.parent.mkdir(parents=True)
        original = self.reference()
        shutil.copyfile(original, old)
        shutil.copyfile(str(original) + ".fai", str(old) + ".fai")
        cache = self.folder / "new_cache"
        resources, _ = REF.resource_plan({"assembly": "hg38", "cache_root": str(cache), "ref_genome": str(old)},
                                        self.data, "selected", use_cache_root=True)
        with self.assertRaisesRegex(ValueError, "lacks verified pinned UCSC provenance"):
            REF.preflight(resources)
        self.assertFalse(cache.exists())

    def test_missing_custom_references_and_chain_do_not_fall_back_to_managed_defaults(self):
        cache = self.folder / "new_cache"
        for key in ("ref_genome", "gofcards_hg19_fasta", "gofcards_hg38_fasta", "gofcards_hg19_to_hg38_chain"):
            with self.subTest(key=key):
                config = {"assembly": "hg38", "cache_root": str(cache), key: str(self.folder / "missing_original.fa")}
                resources, _ = REF.resource_plan(config, self.data, use_cache_root=True)
                with self.assertRaisesRegex(ValueError, "Explicit reference resource is missing"):
                    REF.preflight(resources)
                self.assertFalse(cache.exists())
                self.assertFalse(self.data.exists())

    def test_retained_legacy_index_archives_receipts_staging_and_locks_block_redirection(self):
        old = self.data / "reference/ucsc.hg38.fasta"
        old.parent.mkdir(parents=True)
        cache = self.folder / "new_cache"
        config = {"assembly": "hg38", "cache_root": str(cache), "ref_genome": str(old)}
        for suffix in (".fai", ".fai.partial", ".gzi", ".gz", ".gz.partial", ".partial", ".partial.fai",
                       ".partial.fai.partial", ".source.json", ".source.json.partial", ".install.lock"):
            with self.subTest(suffix=suffix):
                retained = Path(str(old) + suffix)
                retained.write_bytes(b"retained legacy state\n")
                with self.assertRaisesRegex(ValueError, "retained index/source/staging state"):
                    REF.resource_plan(config, self.data, "selected", use_cache_root=True)
                self.assertEqual(retained.read_bytes(), b"retained legacy state\n")
                self.assertFalse(cache.exists())
                retained.unlink()

    def test_dangling_old_default_symlink_and_new_orphan_index_remain_fail_closed(self):
        old = self.data / "reference/ucsc.hg38.fasta"
        old.parent.mkdir(parents=True)
        old.symlink_to(self.folder / "missing_target.fa")
        cache = self.folder / "new_cache"
        resources, _ = REF.resource_plan({"assembly": "hg38", "cache_root": str(cache), "ref_genome": str(old)},
                                        self.data, "selected", use_cache_root=True)
        with self.assertRaisesRegex(ValueError, "Missing, empty, or non-file"):
            REF.preflight(resources)
        self.assertTrue(old.is_symlink())
        new = cache / "reference/ucsc.hg38.fasta"
        new.parent.mkdir(parents=True)
        index = Path(str(new) + ".fai")
        index.write_text("retained orphan index\n")
        resources, _ = REF.resource_plan({"assembly": "hg38", "cache_root": str(cache)}, self.data,
                                        "selected", use_cache_root=True)
        with self.assertRaisesRegex(ValueError, "Orphan FASTA index"):
            REF.preflight(resources)
        self.assertEqual(index.read_text(), "retained orphan index\n")
        self.assertFalse(new.exists())

    def test_forward_chain_routing_keeps_explicit_direction_checks(self):
        chain = self.folder / "selected_reverse.chain"
        original = "chain 1 chr1 14 + 0 1 chr1 12 + 0 1 1\n1\n"
        chain.write_text(original)
        cache = self.folder / "new_cache"
        resources, _ = REF.resource_plan({"assembly": "GRCh37", "cache_root": str(cache),
                                         "gofcards_hg19_to_hg38_chain": str(chain)}, self.data, use_cache_root=True)
        with self.assertRaisesRegex(ValueError, "Expected the forward UCSC"):
            REF.preflight(resources)
        self.assertEqual(chain.read_text(), original)
        self.assertFalse(cache.exists())

    def test_failed_cache_root_install_keeps_config_until_all_resources_succeed(self):
        sources, bodies, archives = self.sources()
        cache = self.folder / "new_cache"
        config_path = self.folder / "config.yaml"
        original = f"assembly: hg38\ncache_root: {cache}\nref_genome: {self.data}/reference/ucsc.hg38.fasta\n# preserve\n"
        config_path.write_text(original)
        argv = ["reference_resources.py", "install", str(config_path), "--data-dir", str(self.data),
                "--scope", "gofcards", "--write-config", "--use-cache-root"]

        def interrupted(request, **_):
            if request.full_url == sources["hg38"]["url"]:
                return Response(archives[request.full_url])
            raise OSError("fixture interruption after first resource")

        with mock.patch.dict(REF.SOURCES, sources), mock.patch.object(sys, "argv", argv), \
                contextlib.redirect_stderr(io.StringIO()):
            with mock.patch.object(REF, "urlopen", side_effect=interrupted):
                self.assertEqual(REF.main(), 1)
            self.assertEqual(config_path.read_text(), original)
            selected = cache / "reference/ucsc.hg38.fasta"
            self.assertEqual(selected.read_bytes(), bodies["hg38"])
            self.assertTrue(REF.valid_index(selected, "hg38"))
            self.assertFalse(self.data.exists())
            with mock.patch.object(REF, "urlopen", side_effect=lambda request, **_: Response(archives[request.full_url])):
                self.assertEqual(REF.main(), 0)
        updated = config_path.read_text()
        self.assertIn("# preserve", updated)
        self.assertIn(str(cache / "reference/ucsc.hg19.fasta"), updated)
        self.assertIn(str(cache / "liftOver/hg19ToHg38.over.chain"), updated)
        self.assertFalse(self.data.exists())

    def test_prerequisite_and_install_planners_agree_on_cache_root_destinations(self):
        cache = self.folder / "new_cache"
        vep = cache / "VEP_cache"
        for assembly in ("GRCh37", "GRCh38"):
            for relative in ("info.txt", "chr_synonyms.txt", "1/all_vars.gz", "1/all_vars.gz.csi"):
                path = vep / "homo_sapiens_merged" / f"116_{assembly}" / relative
                path.parent.mkdir(parents=True, exist_ok=True)
                path.write_text(f"assembly\t{assembly}\nspecies\thomo_sapiens\nvar_type\ttabix\n"
                                if relative == "info.txt" else "fixture\n")
        config = {"assembly": "GRCh37", "cache_root": str(cache), "vep_cache_dir": str(vep),
                  "vep_installed_version": "116", "ref_genome": str(self.data / "reference/ucsc.hg19.fasta")}
        expected, _ = REF.resource_plan(config, self.data, use_cache_root=True)
        with mock.patch.object(PRE.REFERENCES, "preflight", wraps=PRE.REFERENCES.preflight) as preflight:
            PRE.validate_gofcards(config, main_preflight=True, allow_pinned_reference_defaults=True,
                                 source_root=self.data.parent)
        self.assertEqual([row["path"] for row in preflight.call_args.args[0]], [row["path"] for row in expected])
        self.assertFalse((cache / "reference").exists())
        self.assertFalse((cache / "liftOver").exists())
        self.assertFalse(self.data.exists())

    def test_invalid_cache_roots_are_rejected_only_when_installer_mode_opts_in(self):
        for value in ("relative_cache", True, "/", "/paedyl01/original/cache", str(self.folder / "cache with spaces")):
            with self.subTest(value=value):
                config = {"assembly": "hg38", "cache_root": value}
                with self.assertRaises(ValueError):
                    REF.resource_plan(config, self.data, "selected", use_cache_root=True)
                REF.preflight(REF.resource_plan(config, self.data, "selected")[0])
                self.assertFalse(self.data.exists())


    def test_missing_original_dragen_is_never_a_download_destination(self):
        path = self.folder / "DRAGEN/original_reference.fa"
        resources, _ = REF.resource_plan({"assembly": "hg38", "ref_genome": str(path)}, self.data)
        with mock.patch.object(REF, "urlopen", side_effect=AssertionError("Unexpected download")):
            with self.assertRaisesRegex(ValueError, "restore the selected original input"):
                REF.preflight(resources)
        self.assertFalse(path.parent.exists())
        self.assertFalse(self.data.exists())

    def test_supplied_selected_reference_is_reused_under_its_original_name(self):
        path = self.reference(name="DRAGEN.original_reference.fa")
        before = path.read_bytes()
        resources, updates = REF.resource_plan({"assembly": "hg38", "ref_genome": str(path)}, self.data, "both")
        self.assertEqual(updates["gofcards_hg38_fasta"], str(path))
        selected = next(resource for resource in resources if resource["path"] == path)
        with mock.patch.object(REF, "urlopen", side_effect=AssertionError("Unexpected download")):
            REF.preflight([selected])
            REF.install_resource(selected)
        self.assertEqual(path.read_bytes(), before)
        self.assertFalse(REF.receipt_path(path).exists())

    def test_wrong_supplied_assembly_fails_before_writes_or_downloads(self):
        path = self.reference("hg19", "misleading.hg38.fasta")
        before = path.read_bytes()
        resources, _ = REF.resource_plan({"assembly": "hg19", "gofcards_hg38_fasta": str(path)}, self.data)
        with mock.patch.object(REF, "urlopen", side_effect=AssertionError("Unexpected download")):
            with self.assertRaisesRegex(ValueError, "required hg38 chromosome 1"):
                REF.preflight(resources)
        self.assertEqual(path.read_bytes(), before)
        self.assertFalse(self.data.exists())

    def test_same_file_cannot_represent_both_assemblies(self):
        path = self.reference()
        alias = self.folder / "alias.fa"
        alias.symlink_to(path.name)
        with self.assertRaisesRegex(ValueError, "distinct assembly-correct"):
            REF.resource_plan({"assembly": "hg38", "gofcards_hg19_fasta": str(alias),
                               "gofcards_hg38_fasta": str(path)}, self.data)

    def test_unproven_existing_default_is_retained(self):
        path = self.data / "reference/ucsc.hg38.fasta"
        path.parent.mkdir(parents=True)
        path.write_text(">chr1\nAAAAAAAAAAAAAA\n")
        before = path.read_bytes()
        resources, _ = REF.resource_plan({"assembly": "hg38"}, self.data, "selected")
        with self.assertRaisesRegex(ValueError, "lacks verified pinned UCSC provenance"):
            REF.preflight(resources)
        self.assertEqual(path.read_bytes(), before)
        self.assertEqual(list(path.parent.iterdir()), [path])

    def test_both_real_indexes_chain_and_source_receipts_install_and_reuse(self):
        resources, updates, bodies = self.installed()
        for resource in resources:
            path = resource["path"]
            self.assertEqual(path.read_bytes(), bodies[resource["assembly"]])
            self.assertTrue(Path(str(path) + ".gz").is_file())
            receipt = json.loads(REF.receipt_path(path).read_text())
            self.assertEqual(receipt["source"], REF.SOURCES[resource["assembly"]])
            self.assertEqual(receipt["sha256"], hashlib.sha256(path.read_bytes()).hexdigest())
            if resource["kind"] == "fasta":
                self.assertTrue(REF.valid_index(path, resource["assembly"]))
        self.assertEqual(updates["gofcards_hg19_to_hg38_chain"], str(self.data / "liftOver/hg19ToHg38.over.chain"))
        with mock.patch.object(REF, "urlopen", side_effect=AssertionError("Unexpected download")):
            REF.preflight(resources)
            for resource in resources:
                REF.install_resource(resource)

    def test_modified_managed_reference_is_not_overwritten(self):
        resources, _, _ = self.installed()
        path = resources[0]["path"]
        changed = path.read_bytes().replace(b"A", b"C", 1)
        path.write_bytes(changed)
        with self.assertRaisesRegex(ValueError, "checksum mismatch"):
            REF.preflight(resources)
        self.assertEqual(path.read_bytes(), changed)

    def test_lost_receipt_recovers_only_from_the_exact_pinned_archive(self):
        resources, _, _ = self.installed()
        resource = resources[0]
        REF.receipt_path(resource["path"]).unlink()
        with mock.patch.object(REF, "urlopen", side_effect=AssertionError("Unexpected download")):
            REF.preflight(resources)
            REF.install_resource(resource)
        self.assertTrue(REF.receipt_path(resource["path"]).is_file())

    def test_receipt_for_a_different_source_is_rejected(self):
        resources, _, _ = self.installed()
        receipt = REF.receipt_path(resources[0]["path"])
        metadata = json.loads(receipt.read_text())
        metadata["source"]["url"] += ".another-release"
        receipt.write_text(json.dumps(metadata))
        with self.assertRaisesRegex(ValueError, "source identity mismatch"):
            REF.preflight(resources)

    def test_missing_or_wrong_index_is_repaired_without_reference_changes(self):
        for index in (False, True):
            with self.subTest(index=index):
                path = self.reference(name=f"original_{index}.fa", index=index)
                if index:
                    Path(str(path) + ".fai").write_text("chr1\t12\t6\t12\t13\n")
                before = path.read_bytes()
                resource = {"path": path, "assembly": "hg38", "kind": "fasta", "managed": False}
                REF.preflight([resource])
                REF.install_resource(resource)
                self.assertTrue(REF.valid_index(path, "hg38"))
                self.assertEqual(path.read_bytes(), before)

    def test_forward_chain_rejects_reverse_direction(self):
        path = self.folder / "chain"
        path.write_text("chain 1 chr1 14 + 0 1 chr1 12 + 0 1 1\n1\n")
        with self.assertRaisesRegex(ValueError, "forward UCSC"):
            REF.validate_chain(path)

    def test_existing_wrong_archive_is_never_replaced(self):
        sources, _, _ = self.sources()
        path = self.folder / "reference.gz"
        path.write_bytes(b"existing wrong source")
        with mock.patch.object(REF, "urlopen", side_effect=AssertionError("Unexpected download")):
            with self.assertRaisesRegex(ValueError, "Existing UCSC source archive checksum mismatch"):
                REF.download(sources["hg19"], path)
        self.assertEqual(path.read_bytes(), b"existing wrong source")

    def test_interrupted_download_resumes_at_exact_partial_length(self):
        sources, _, archives = self.sources()
        source = sources["hg19"]
        body = archives[source["url"]]
        path = self.folder / "reference.gz"
        with mock.patch.object(REF, "urlopen", return_value=InterruptedResponse(body)):
            with self.assertRaisesRegex(OSError, "interrupted"):
                REF.download(source, path)
        partial = Path(str(path) + ".partial")
        self.assertEqual(partial.read_bytes(), body[:8])
        self.assertFalse(path.exists())

        def resume(request, **_):
            self.assertEqual(request.get_header("Range"), "bytes=8-")
            return Response(body[8:], 206, {"Content-Range": f"bytes 8-{len(body)-1}/{len(body)}"})

        with mock.patch.object(REF, "urlopen", side_effect=resume):
            REF.download(source, path)
        self.assertEqual(path.read_bytes(), body)
        self.assertFalse(partial.exists())

    def test_ignored_range_restarts_the_partial(self):
        sources, _, archives = self.sources()
        source = sources["hg19"]
        path = self.folder / "reference.gz"
        Path(str(path) + ".partial").write_bytes(b"truncated")
        body = archives[source["url"]]
        with mock.patch.object(REF, "urlopen", return_value=Response(body)):
            REF.download(source, path)
        self.assertEqual(path.read_bytes(), body)

    def test_wrong_range_does_not_append_or_publish(self):
        sources, _, archives = self.sources()
        path = self.folder / "reference.gz"
        partial = Path(str(path) + ".partial")
        partial.write_bytes(b"partial")
        with mock.patch.object(REF, "urlopen", return_value=Response(b"wrong", 206, {"Content-Range": "bytes 0-4/5"})):
            with self.assertRaisesRegex(ValueError, "Incorrect UCSC download Content-Range"):
                REF.download(sources["hg19"], path)
        self.assertEqual(partial.read_bytes(), b"partial")
        self.assertFalse(path.exists())

    def test_verified_complete_partial_needs_no_download(self):
        sources, _, archives = self.sources()
        source = sources["hg19"]
        path = self.folder / "reference.gz"
        Path(str(path) + ".partial").write_bytes(archives[source["url"]])
        with mock.patch.object(REF, "urlopen", side_effect=AssertionError("Unexpected download")):
            REF.download(source, path)
        self.assertEqual(path.read_bytes(), archives[source["url"]])

    def test_corrupt_download_is_retained_as_partial_only(self):
        sources, _, archives = self.sources()
        source = sources["hg19"]
        path = self.folder / "reference.gz"
        corrupted = b"x" * source["bytes"]
        with mock.patch.object(REF, "urlopen", return_value=Response(corrupted)):
            with self.assertRaisesRegex(ValueError, "checksum/size mismatch"):
                REF.download(source, path)
        self.assertFalse(path.exists())
        self.assertEqual(Path(str(path) + ".partial").read_bytes(), corrupted)

    def test_symlink_outside_policy_is_rejected_before_target_inspection(self):
        path = self.folder / "outside"
        path.symlink_to("/not-permitted/private-reference")
        with self.assertRaisesRegex(ValueError, "symlink leaves permitted inspection roots"):
            REF.checked_path(str(path))

    def test_main_preflight_defers_only_pinned_defaults_and_still_checks_opposite_vep(self):
        cache = self.folder / "cache"
        for assembly in ("GRCh37", "GRCh38"):
            for relative in ("info.txt", "chr_synonyms.txt", "1/all_vars.gz", "1/all_vars.gz.csi"):
                path = cache / "homo_sapiens_merged" / f"113_{assembly}" / relative
                path.parent.mkdir(parents=True, exist_ok=True)
                path.write_text(f"species\thomo_sapiens\nassembly\t{assembly}\nvar_type\ttabix\n"
                                if relative == "info.txt" else "fixture\n")
        config = {"assembly": "hg19", "vep_cache_dir": str(cache), "vep_installed_version": "113"}
        with mock.patch.dict(os.environ, {"PRIVA_SMOKE_SKIP_VEP": "0"}):
            PRE.validate_gofcards(config, main_preflight=True, allow_pinned_reference_defaults=True,
                                 source_root=self.data.parent)
            PRE.validate_gofcards(config, allow_pinned_reference_defaults=True, source_root=self.data.parent)
            self.assertFalse(self.data.exists())
            config["gofcards_hg38_fasta"] = str(self.folder / "missing_original.fa")
            with self.assertRaisesRegex(ValueError, "Explicit reference resource is missing"):
                PRE.validate_gofcards(config, main_preflight=True, allow_pinned_reference_defaults=True,
                                     source_root=self.data.parent)
            config.pop("gofcards_hg38_fasta")
            (cache / "homo_sapiens_merged/113_GRCh38/1/all_vars.gz").unlink()
            with self.assertRaisesRegex(ValueError, "GoFCards GRCh38 merged VEP"):
                PRE.validate_gofcards(config, main_preflight=True, allow_pinned_reference_defaults=True,
                                     source_root=self.data.parent)
            with self.assertRaisesRegex(ValueError, "GoFCards GRCh38 merged VEP"):
                PRE.validate_gofcards(config, allow_pinned_reference_defaults=True, source_root=self.data.parent)


    def test_configuration_paths_are_quoted_and_unrelated_comments_preserved(self):
        original = 'assembly: hg38\nref_genome: null  # selected reference\nthreads: 7  # keep\n'
        path = str(self.folder / "reference & special # label.fa")
        updated = REF.render_config_values(original, {"ref_genome": path, "gofcards_hg38_fasta": path})
        self.assertIn("# selected reference", updated)
        self.assertIn("threads: 7  # keep", updated)
        config = REF.yaml.safe_load(updated)
        self.assertEqual(config["ref_genome"], path)
        self.assertEqual(config["gofcards_hg38_fasta"], path)
        self.assertEqual(REF.render_config_values(updated, {"ref_genome": path, "gofcards_hg38_fasta": path}), updated)

    def test_json_configuration_and_concurrent_configuration_changes(self):
        original = '{"assembly": "hg38", "threads": 7}\n'
        updated = REF.render_config_values(original, {"ref_genome": str(self.folder / "ref.fa")})
        self.assertEqual(json.loads(updated)["threads"], 7)
        path = self.folder / "config.json"
        path.write_text(original)
        REF.write_config_values(path, original, updated)
        self.assertEqual(path.read_text(), updated)
        changed = '{"assembly": "hg19", "threads": 9}\n'
        path.write_text(changed)
        with self.assertRaisesRegex(ValueError, "Configuration changed"):
            REF.write_config_values(path, original, updated)
        self.assertEqual(path.read_text(), changed)

    def test_download_failure_never_writes_configured_paths(self):
        config = self.folder / "config.yaml"
        original = "assembly: hg38\nref_genome: null  # keep until successful\n"
        config.write_text(original)
        argv = ["reference_resources.py", "install", str(config), "--data-dir", str(self.data), "--write-config"]
        with mock.patch.object(sys, "argv", argv), \
             mock.patch.object(REF, "urlopen", side_effect=OSError("fixture network failure")), \
             contextlib.redirect_stderr(io.StringIO()):
            self.assertEqual(REF.main(), 1)
        self.assertEqual(config.read_text(), original)


if __name__ == "__main__":
    unittest.main()
