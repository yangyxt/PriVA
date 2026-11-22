"""Small synthetic download/retry and explicit source identity regressions."""
import argparse
from dataclasses import replace
import importlib
import json
import os
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest
from unittest import mock

import fake_priva_binary

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
INPUTS = importlib.import_module("mechanism_source_input")
BUILDERS = [importlib.import_module(name) for name in (
    "build_gene_pathogenic_mechanism_cache", "build_gene_nonlof_mechanism_cache")]
EREPO = "ClinVar Variation Id\tApplied Evidence Codes (Met)\n123\tPS3\n"
DOSAGE = ("# published metadata\n#Gene Symbol\tGenomic Location\tHaploinsufficiency Score\t"
          "Haploinsufficiency Description\tHaploinsufficiency Disease ID\tTriplosensitivity Score\n"
          "GENE\tchr1:1-2\t3\tSufficient evidence\tMONDO:0000001\t0\n")
G2P = "gene symbol,disease name,molecular mechanism\nGENE,Condition,gain of function\n"


class ClinGenRecovery(unittest.TestCase):
    def setUp(self):
        temp = tempfile.TemporaryDirectory(prefix="clingen_retry_", dir=os.environ.get("PRIVA_TEST_TMPDIR", os.environ.get("TMPDIR")))
        self.addCleanup(temp.cleanup)
        self.folder = Path(temp.name)
        self.cache = self.folder / "cache"
        self.data = self.folder / "bundled"
        self.cache.mkdir()
        self.erepo = self.cache / "clingen_erepo.tsv"
        self.dosage = self.cache / "dosage.hg38.tsv"
        self.map = self.cache / "clingen_map.hg38.pkl.gz"
        self.config = self.folder / "config.json"
        self.config.write_text(json.dumps({"assembly": "hg38", "cache_root": str(self.cache),
            "clinvar_vcf": str(self.folder / "clinvar.vcf"),
            "clingen_map": str(self.map), "gene_dosage_sensitivity": str(self.dosage),
            "clingen_dosage_grch38_url": "https://selected.invalid/GRCh38.tsv"}))
        (self.folder / "erepo.tsv").write_text(EREPO)
        (self.folder / "dosage.tsv").write_text(DOSAGE)
        self.env = {key: value for key, value in os.environ.items() if not key.startswith("PRIVA_FORCE_")}
        self.env.update(TEST_PYTHON=sys.executable, TEST_FOLDER=str(self.folder), TEST_CACHE=str(self.cache),
                        TEST_DATA=str(self.data), TEST_CONFIG=str(self.config), TEST_MODE="success",
                        **fake_priva_binary.environment(self.folder))


    def test_download_validator_rejects_incomplete_source_rows(self):
        for kind, payload in (("clingen_erepo", EREPO), ("clingen_dosage", DOSAGE)):
            with self.subTest(kind=kind):
                target = self.folder / "incomplete.tsv"
                target.write_text(payload + "truncated\n")
                with self.assertRaisesRegex(ValueError, "incomplete"):
                    INPUTS.validate_source_payload(target, kind)

    def test_shell_download_validator_requires_runtime_exact_headers(self):
        path = self.folder / "aliases.tsv"
        for kind, payload in (("clingen_erepo", "clinvar_variation_id\tapplied_evidence_codes_met\n123\tPS3\n"),
                              ("clingen_dosage_runtime", "symbol\thi_score\tts_score\nGENE\t3\t0\n")):
            path.write_text(payload)
            with self.assertRaisesRegex(ValueError, "exact columns"):
                INPUTS.validate_source_payload(path, kind)
        # Flexible builder inputs remain distinct from strict shell-runtime inputs.
        INPUTS.validate_source_payload(path, "clingen_dosage")
        path.write_text(EREPO.replace("123\t", "not-an-id\t"))
        with self.assertRaisesRegex(ValueError, "nonnumeric"):
            INPUTS.validate_source_payload(path, "clingen_erepo")


class SelectedSourceURLs(unittest.TestCase):
    def setUp(self):
        temp = tempfile.TemporaryDirectory(prefix="source_urls_", dir=os.environ.get("PRIVA_TEST_TMPDIR", os.environ.get("TMPDIR")))
        self.addCleanup(temp.cleanup)
        self.folder = Path(temp.name)
        for module in BUILDERS:
            module.configure_source_urls(argparse.Namespace())
        self.addCleanup(lambda: [m.configure_source_urls(argparse.Namespace()) for m in BUILDERS])

    def fetch(self, module, source, previous):
        return module.fetch_static_source(source=source, raw_dir=self.folder, previous_meta=previous,
            force=False, timeout=1, retries=1, proxy_url="", download_tool="urllib")

    def test_default_and_named_cli_selection_retain_parser_order_and_refresh(self):
        for module in BUILDERS:
            with self.subTest(module=module.__name__):
                default = [(s.name, s.url, s.parser, s.refresh_days) for s in module.SOURCES]
                with mock.patch.object(sys, "argv", [module.__file__, "--g2p-ddg2p-url", "https://selected.invalid/g2p.csv"]):
                    args = module.parse_args()
                module.configure_source_urls(args)
                self.assertEqual(module.SOURCES[0].url, "https://selected.invalid/g2p.csv")
                self.assertEqual([(s.name, s.parser, s.refresh_days) for s in module.SOURCES], [(n, p, d) for n, _, p, d in default])
                module.configure_source_urls(argparse.Namespace())
                self.assertEqual([(s.name, s.url, s.parser, s.refresh_days) for s in module.SOURCES], default)

    def test_changed_url_fetches_even_when_old_raw_is_fresh_and_records_identity(self):
        for module in BUILDERS:
            with self.subTest(module=module.__name__):
                old = module.SOURCES[0]
                source = replace(old, url="https://selected.invalid/new.csv")
                raw = self.folder / source.raw_filename
                raw.write_text("old source bytes\n")
                previous = {source.name: {"url": old.url, "raw_path": str(raw), "sha256": module.sha256_file(raw)}}
                def download(url, output, **kwargs):
                    self.assertEqual(url, source.url)
                    output.write_text(G2P)
                with mock.patch.object(module, "download_url", side_effect=download) as fetch:
                    metadata = self.fetch(module, source, previous)
                self.assertEqual(fetch.call_count, 1)
                self.assertEqual(metadata["url"], source.url)
                self.assertEqual(raw.read_text(), G2P)
                with mock.patch.object(module, "download_url", side_effect=AssertionError("unexpected download")):
                    reused = self.fetch(module, source, {source.name: metadata})
                self.assertEqual(reused["status"], "skipped_fresh")

    def test_failed_changed_url_or_incompatible_schema_never_relabels_old_raw(self):
        for module in BUILDERS:
            old = module.SOURCES[0]
            source = replace(old, url="https://selected.invalid/new.csv")
            raw = self.folder / source.raw_filename
            for mode in ("network", "schema", "unreceipted"):
                with self.subTest(module=module.__name__, mode=mode):
                    raw.write_text("previous bytes\n")
                    previous = {} if mode == "unreceipted" else {source.name: {"url": old.url, "raw_path": str(raw)}}
                    def download(url, output, **kwargs):
                        if mode in ("network", "unreceipted"):
                            output.write_text("partial\n")
                            raise OSError("injected transfer failure")
                        output.write_text("wrong,columns\n1,2\n")
                    with mock.patch.object(module, "download_url", side_effect=download):
                        with self.assertRaisesRegex(RuntimeError, "prior bytes preserved"):
                            self.fetch(module, source, previous)
                    self.assertEqual(raw.read_text(), "previous bytes\n")

    def test_unchanged_legacy_default_fresh_bytes_keep_existing_reuse(self):
        for module in BUILDERS:
            source = module.SOURCES[0]
            raw = self.folder / source.raw_filename
            raw.write_text("opaque existing default source\n")
            with mock.patch.object(module, "download_url", side_effect=AssertionError("unexpected download")):
                metadata = self.fetch(module, source, {})
            self.assertEqual(metadata["status"], "skipped_fresh")

    def test_old_manifest_cannot_relabel_raw_replaced_by_failed_other_source_build(self):
        for module in BUILDERS:
            source = module.SOURCES[0]
            raw = self.folder / source.raw_filename
            raw.write_text(G2P)
            previous = {source.name: {"url": source.url, "raw_path": str(raw), "sha256": module.sha256_file(raw)}}
            raw.write_text("new B bytes from interrupted changed-source build\n")
            with mock.patch.object(module, "download_url", side_effect=OSError("injected unavailable A")) as fetch:
                with self.assertRaisesRegex(RuntimeError, "prior bytes preserved"):
                    self.fetch(module, source, previous)
            self.assertEqual(fetch.call_count, 1)
            self.assertEqual(raw.read_text(), "new B bytes from interrupted changed-source build\n")
            with mock.patch.object(module, "download_url", side_effect=lambda url, output, **kwargs: output.write_text(G2P)):
                metadata = self.fetch(module, source, previous)
            self.assertEqual(metadata["url"], source.url)
            self.assertEqual(raw.read_text(), G2P)

    def test_identical_original_panelapp_download_restores_replaced_raw_bytes(self):
        for module in BUILDERS:
            source = module.PANELAPP_SPEC
            raw = self.folder / source.raw_filename
            responses = [{"results": [{"id": 7}], "next": None}, {"id": 7, "genes": []}]
            with mock.patch.object(module, "iso_now", return_value="2026-10-18T00:00:00+00:00"), \
                 mock.patch.object(module, "fetch_json", side_effect=responses):
                original = module.fetch_panelapp(self.folder, {}, True, 1, 1, "", "urllib", None)
            payload = raw.read_bytes()
            raw.write_bytes(b"other source bytes from failed build\n")
            with mock.patch.object(module, "iso_now", return_value="2026-10-18T00:00:00+00:00"), \
                 mock.patch.object(module, "fetch_json", side_effect=responses):
                restored = module.fetch_panelapp(self.folder, {source.name: original}, False, 1, 1, "", "urllib", None)
            self.assertEqual(raw.read_bytes(), payload)
            self.assertEqual(restored["url"], source.url)

    def test_panelapp_override_uses_selected_list_and_detail_host(self):
        for module in BUILDERS:
            with self.subTest(module=module.__name__):
                module.configure_source_urls(argparse.Namespace(panelapp_panels_url="https://selected.invalid/api/panels/"))
                # Only the non-LOF builder exposes PanelApp as an active input.
                if not module.__name__.endswith("nonlof_mechanism_cache"):
                    self.assertEqual(module.PANELAPP_SPEC.url, module.DEFAULT_PANELAPP_SPEC.url)
                    continue
                source = module.PANELAPP_SPEC
                raw = self.folder / source.raw_filename
                raw.write_text("previous PanelApp\n")
                previous = {source.name: {"url": module.DEFAULT_PANELAPP_SPEC.url, "raw_path": str(raw)}}
                calls = []
                def fetch(url, **kwargs):
                    calls.append(url)
                    return {"results": [{"id": 7}], "next": None} if len(calls) == 1 else {"id": 7, "genes": []}
                with mock.patch.object(module, "fetch_json", side_effect=fetch):
                    metadata = module.fetch_panelapp(self.folder, previous, False, 1, 1, "", "urllib", None)
                self.assertEqual(calls, [source.url, "https://selected.invalid/api/panels/7/?format=json"])
                self.assertEqual(metadata["url"], source.url)

    def test_bad_panelapp_schema_or_detail_failure_preserves_previous_file(self):
        module = BUILDERS[1]
        module.configure_source_urls(argparse.Namespace(panelapp_panels_url="https://selected.invalid/panels/"))
        raw = self.folder / module.PANELAPP_SPEC.raw_filename
        for mode in ("schema", "detail", "nested"):
            raw.write_text("previous PanelApp\n")
            responses = [{"wrong": []}] if mode == "schema" else [{"results": [{"id": 7}], "next": None},
                {"id": 7, "genes": [None]} if mode == "nested" else OSError("detail failed")]
            with mock.patch.object(module, "fetch_json", side_effect=responses):
                with self.assertRaisesRegex(RuntimeError, "prior bytes preserved"):
                    module.fetch_panelapp(self.folder, {}, False, 1, 1, "", "urllib", None)
            self.assertEqual(raw.read_text(), "previous PanelApp\n")

    def test_panelapp_return_to_default_url_keeps_strict_transition_failure(self):
        module = BUILDERS[1]
        source = module.PANELAPP_SPEC
        raw = self.folder / source.raw_filename
        raw.write_text("previous custom PanelApp\n")
        previous = {source.name: {"url": "https://previous.invalid/panels/", "raw_path": str(raw)}}
        with mock.patch.object(module, "fetch_json", side_effect=[{"results": [{"id": 7}], "next": None}, OSError("detail failed")]):
            with self.assertRaisesRegex(RuntimeError, "prior bytes preserved"):
                module.fetch_panelapp(self.folder, previous, False, 1, 1, "", "urllib", None)
        self.assertEqual(raw.read_text(), "previous custom PanelApp\n")

    def test_clinvar_only_refresh_cannot_relabel_changed_source(self):
        module = BUILDERS[1]
        old = module.SOURCES[0]
        selected = replace(old, url="https://selected.invalid/new.csv")
        raw = self.folder / old.raw_filename
        raw.write_text(G2P)
        with self.assertRaisesRegex(ValueError, "cannot change selected source"):
            module.carry_forward_source_metadata(selected, self.folder, {old.name: {"url": old.url}})

    def test_orphadata_override_requires_complete_original_xml_schema(self):
        module = BUILDERS[0]
        module.configure_source_urls(argparse.Namespace(orphadata_gene_disease_url="https://selected.invalid/product6.xml"))
        source = next(s for s in module.CONDITION_EVIDENCE_SPECS if s.parser == "orphadata")
        raw = self.folder / source.raw_filename
        raw.parent.mkdir(parents=True)
        original = b"previous Product 6 source\n"
        raw.write_bytes(original)
        for payload in ("<html/>", "<JDBOR><DisorderList>", '<JDBOR xmlns="urn:unsupported"><DisorderList/></JDBOR>'):
            with mock.patch.object(module, "download_url", side_effect=lambda url, output, **kwargs: output.write_text(payload)):
                with self.assertRaisesRegex(RuntimeError, "prior bytes preserved"):
                    self.fetch(module, source, {})
            self.assertEqual(raw.read_bytes(), original)
        with mock.patch.object(module, "download_url", side_effect=lambda url, output, **kwargs: output.write_text("<JDBOR><DisorderList/></JDBOR>\n")):
            metadata = self.fetch(module, source, {})
        self.assertEqual(metadata["url"], source.url)


if __name__ == "__main__":
    unittest.main()
