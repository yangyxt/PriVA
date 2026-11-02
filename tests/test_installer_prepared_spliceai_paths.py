"""Prepared SpliceAI2 selection survives cache-root planning and registration."""
import argparse
import copy
import gzip
import hashlib
import importlib.util
import json
import os
from pathlib import Path
import tempfile
import unittest
from unittest import mock

import yaml


ROOT = Path(__file__).resolve().parents[1]
SPEC = importlib.util.spec_from_file_location("prepared_spliceai_paths", ROOT / "scripts/installer_config_paths.py")
PATHS = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(PATHS)


class PreparedSpliceAIPaths(unittest.TestCase):
    def setUp(self):
        temporary = tempfile.TemporaryDirectory(prefix="priva_prepared_paths_",
            dir=os.environ.get("PRIVA_TEST_TMPDIR", os.environ.get("TMPDIR")))
        self.addCleanup(temporary.cleanup)
        self.folder = Path(temporary.name)
        self.root = self.folder / "cache"
        self.root.mkdir()
        self.config = self.folder / "config.yaml"
        self.settings = dict(cache_root=str(self.root), tmp_dir=str(self.folder / "tmp"),
                             assembly="hg38", vep_installed_version="116")
        self.cache = PATHS.spliceai_inspector()
        self.canonical = self.root / "VEP_plugin_cache/SpliceAI2/GRCh38"

    def prepared_pair(self, folder, *, label="selected"):
        """Small real files/receipts; synthetic indexes exercise metadata only."""
        folder.mkdir(parents=True)
        outputs = {}
        for role in ("snv", "indel"):
            path = folder / f"SpliceAI2_GRCh38.{role}.vep.vcf.gz"
            text = ("##fileformat=VCFv4.2\n##assembly=GRCh38\n##SpliceAI2Version=2\n"
                    f"##fixture={label}\n"
                    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
            path.write_bytes(gzip.compress(text.encode(), mtime=0))
            index = Path(str(path) + ".tbi")
            index.write_bytes(b"synthetic index; not an indexing acceptance fixture")
            outputs[role] = dict(path=str(path), index_path=str(index), bytes=path.stat().st_size,
                sha256=hashlib.sha256(path.read_bytes()).hexdigest(),
                index_sha256=hashlib.sha256(index.read_bytes()).hexdigest())
        context = folder / "SpliceAI2_GRCh38.gene_context.tsv"
        context.write_text("chrom\tstart\tend\tstrand\tgene\nchr1\t1\t100\t+\tGENE\n")
        receipt = dict(schema_version=1, resource="SpliceAI2", version="2.0", assembly="GRCh38",
            status="prepared_annotation_only", coverage="complete_published", scientific_qualified=False,
            source_revision=self.cache.REVISION, source_manifest_sha256="a" * 64, outputs=outputs,
            gene_context=dict(path=str(context), sha256=hashlib.sha256(context.read_bytes()).hexdigest()))
        (folder / "adapter.ready.json").write_text(json.dumps(receipt))
        (folder / ".install.lock").touch()
        return {"spliceai_" + role + "_prescore": item["path"] for role, item in outputs.items()}

    def select(self, selected):
        self.settings.update(selected)
        self.config.write_text(yaml.safe_dump(self.settings))

    def plan(self):
        updates = PATHS.portable_updates(self.settings, ROOT)
        return updates, PATHS.adoption_plan(self.settings, updates, ROOT)

    def resolve(self):
        PATHS.resolve_config(argparse.Namespace(config=self.config, source_root=ROOT,
            managed_updates_json=None, check=False, output=None))

    def register(self):
        config = yaml.safe_load(self.config.read_text())
        return self.cache.install(self.config, self.root / "VEP_plugin_cache",
            config["spliceai_snv_prescore"], config["spliceai_indel_prescore"])

    @staticmethod
    def snapshot(folder):
        return {str(path.relative_to(folder)): (path.read_bytes(), PATHS.cache_identity(path))
                for path in folder.rglob("*") if path.is_file()}

    def assert_preserved(self, selected):
        updates, adoptions = self.plan()
        effective = dict(self.settings, **updates)
        self.assertEqual({key: effective[key] for key in selected}, selected)
        self.assertFalse(any(item["kind"] == "prepared_projection" for item in adoptions))
        before = self.snapshot(self.root)
        with mock.patch.object(self.cache, "adopt_completed_projection", side_effect=AssertionError("unexpected relocation")), \
                mock.patch.object(self.cache, "_sha256_file", side_effect=AssertionError("unexpected payload hash")):
            self.resolve()
            self.assertEqual(self.register()["format"], "vcf_spliceai2")
        configured = yaml.safe_load(self.config.read_text())
        self.assertEqual({key: configured[key] for key in selected}, selected)
        self.assertEqual(self.snapshot(self.root), before)

    def test_managed_prepared_selection_is_preserved(self):
        selected = self.prepared_pair(self.root / "VEP_plugin_cache/SpliceAI2/selected_names")
        self.select(selected)
        self.assert_preserved(selected)

    def test_root_alias_preserves_both_alias_and_physical_selected_paths(self):
        selected = self.prepared_pair(self.root / "SpliceAI2/selected_names")
        alias = self.folder / "cache_alias"
        alias.symlink_to(self.root, target_is_directory=True)
        self.settings["cache_root"] = str(alias)
        for aliased in (False, True):
            with self.subTest(selected_through_alias=aliased):
                paths = {key: str(alias / Path(value).relative_to(self.root)) if aliased else value
                         for key, value in selected.items()}
                self.select(paths)
                self.assert_preserved(paths)

    def test_different_canonical_pair_does_not_replace_managed_selection(self):
        self.prepared_pair(self.canonical, label="other release")
        selected = self.prepared_pair(self.root / "VEP_plugin_cache/SpliceAI2/selected_names")
        self.select(selected)
        self.assert_preserved(selected)

    def test_external_pair_retains_receipt_aware_physical_adoption(self):
        source = self.folder / "external"
        selected = self.prepared_pair(source)
        before = self.snapshot(source)
        self.select(selected)
        updates, adoptions = self.plan()
        self.assertEqual(adoptions, [dict(key="spliceai2_prepared_pair", source=str(source),
            destination=str(self.canonical), kind="prepared_projection")])
        for key, value in selected.items():
            self.assertEqual(updates[key], str(self.canonical / Path(value).name))
        self.resolve()
        self.assertEqual(self.register()["format"], "vcf_spliceai2")
        self.assertEqual(self.snapshot(source), before)
        self.assertEqual((self.canonical / "adapter.original.ready.json").read_bytes(),
                         (source / "adapter.ready.json").read_bytes())
        for value in selected.values():
            target = self.canonical / Path(value).name
            self.assertEqual(target.read_bytes(), Path(value).read_bytes())
            self.assertNotEqual(target.stat().st_ino, Path(value).stat().st_ino)

    def test_external_collision_fails_before_config_selection_or_replacement(self):
        self.prepared_pair(self.canonical, label="other release")
        source = self.folder / "external"
        self.select(self.prepared_pair(source))
        config_before = self.config.read_bytes()
        snapshots = {folder: self.snapshot(folder) for folder in (source, self.canonical)}
        with self.assertRaisesRegex(ValueError, "different original receipt"):
            self.resolve()
        self.assertEqual(self.config.read_bytes(), config_before)
        for folder, before in snapshots.items():
            self.assertEqual(self.snapshot(folder), before)

    def test_link_beneath_root_to_external_pair_is_still_external(self):
        selected = self.prepared_pair(self.folder / "external")
        link = self.root / "external_link"
        link.symlink_to(self.folder / "external", target_is_directory=True)
        self.select({key: str(link / Path(value).name) for key, value in selected.items()})
        _, adoptions = self.plan()
        self.assertEqual(adoptions, [dict(key="spliceai2_prepared_pair", source=str(self.folder / "external"),
            destination=str(self.canonical), kind="prepared_projection")])

    def test_missing_adapter_receipt_is_rejected_even_inside_root(self):
        folder = self.root / "selected"
        self.select(self.prepared_pair(folder))
        (folder / "adapter.ready.json").unlink()
        with self.assertRaisesRegex(ValueError, "adjacent adapter.ready.json"):
            self.plan()

    def test_managed_receipt_schema_and_pair_gates_still_reject_registration(self):
        folder = self.root / "selected"
        self.select(self.prepared_pair(folder))
        self.resolve()
        receipt = folder / "adapter.ready.json"
        original = json.loads(receipt.read_text())
        invalid = [dict(original, **{key: value}) for key, value in (
            ("coverage", "validation_subset"), ("status", "building"), ("schema_version", 2),
            ("assembly", "GRCh37"), ("source_revision", "different"), ("source_manifest_sha256", None))]
        swapped = copy.deepcopy(original)
        swapped["outputs"]["snv"], swapped["outputs"]["indel"] = swapped["outputs"]["indel"], swapped["outputs"]["snv"]
        invalid.append(swapped)
        wrong_size = copy.deepcopy(original)
        wrong_size["outputs"]["snv"]["bytes"] += 1
        invalid.append(wrong_size)
        for number, data in enumerate(invalid):
            with self.subTest(invalid_receipt=number):
                receipt.write_text(json.dumps(data))
                before = self.config.read_bytes()
                with self.assertRaises(ValueError):
                    self.register()
                self.assertEqual(self.config.read_bytes(), before)

    def test_external_incomplete_receipt_fails_before_config_selection(self):
        folder = self.folder / "external"
        self.select(self.prepared_pair(folder))
        receipt = folder / "adapter.ready.json"
        data = json.loads(receipt.read_text())
        data["coverage"] = "validation_subset"
        receipt.write_text(json.dumps(data))
        before = self.config.read_bytes()
        with self.assertRaisesRegex(ValueError, "completed full GRCh38"):
            self.resolve()
        self.assertEqual(self.config.read_bytes(), before)
        self.assertFalse(self.canonical.exists())

    def test_managed_pair_must_still_match_configured_assembly(self):
        self.select(self.prepared_pair(self.root / "selected"))
        self.settings["assembly"] = "hg19"
        with self.assertRaisesRegex(ValueError, "assembly conflicts"):
            self.plan()


if __name__ == "__main__":
    unittest.main()
