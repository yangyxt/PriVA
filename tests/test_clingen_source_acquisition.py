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
        self.config.write_text(json.dumps({"assembly": "hg38", "clinvar_vcf": str(self.folder / "clinvar.vcf"),
            "clingen_map": str(self.map), "gene_dosage_sensitivity": str(self.dosage),
            "clingen_dosage_grch38_url": "https://selected.invalid/GRCh38.tsv"}))
        (self.folder / "erepo.tsv").write_text(EREPO)
        (self.folder / "dosage.tsv").write_text(DOSAGE)
        self.env = {key: value for key, value in os.environ.items() if not key.startswith("PRIVA_FORCE_")}
        self.env.update(TEST_PYTHON=sys.executable, TEST_FOLDER=str(self.folder), TEST_CACHE=str(self.cache),
                        TEST_DATA=str(self.data), TEST_CONFIG=str(self.config), TEST_MODE="success")

    def run_shell(self, body='ClinGen_deploy "$TEST_CONFIG"', **env):
        prelude = r'''
source "$1" >/dev/null || exit
DATA_DIR="$TEST_DATA"
read_yaml() { "$TEST_PYTHON" -c 'import json,sys; v=json.load(open(sys.argv[1])).get(sys.argv[2]); print("null" if v is None else v)' "$1" "$2"; }
update_yaml() { printf 'UPDATE %s %s\n' "$2" "$3"; }
wget() {
    printf 'FETCH %s\n' "$*"
    local target="${@: -1}"
    case $TEST_MODE in
        fail) printf 'plausible but unfinished\n' > "$target"; return 71 ;;
        empty) : > "$target"; return 0 ;;
        invalid) printf '<html>temporary service error</html>\n' > "$target"; return 0 ;;
    esac
    if [[ "$2" == *classifications* ]]; then cp "$TEST_FOLDER/erepo.tsv" "$target"
    else cp "$TEST_FOLDER/dosage.tsv" "$target"; fi
}
python() {
    if [[ $1 == */prepare_clingen_map.py ]]; then
        printf 'MAP %s\n' "$*"
        while [[ $# -gt 0 ]]; do
            if [[ $1 == -o ]]; then printf 'scientific-map-fixture\n' > "$2"; return; fi
            shift
        done
        return 91
    fi
    "$TEST_PYTHON" "$@"
}
'''
        return subprocess.run(["bash", "--noprofile", "--norc", "-c", prelude + body, "test",
                               str(ROOT / "scripts/install_utils.sh")], text=True, capture_output=True,
                              env=dict(self.env, **env), timeout=40)

    def test_interrupted_erepo_retry_acquires_without_canonical_partial(self):
        self.dosage.write_text(DOSAGE)
        INPUTS.source_receipt(self.dosage, "https://selected.invalid/GRCh38.tsv", record=True)
        failed = self.run_shell(TEST_MODE="fail")
        self.assertNotEqual(failed.returncode, 0)
        self.assertFalse(self.erepo.exists())
        self.assertFalse(self.map.exists())
        self.assertTrue(list(self.cache.glob("clingen_erepo.tsv.download.*")))
        retried = self.run_shell()
        self.assertEqual(retried.returncode, 0, retried.stderr)
        self.assertEqual(self.erepo.read_text(), EREPO)
        self.assertEqual(retried.stdout.count("FETCH "), 1)
        self.assertIn("UPDATE clingen_map ", retried.stdout)

    def test_interrupted_dosage_retry_preserves_erepo_and_uses_selected_url(self):
        self.erepo.write_text(EREPO)
        failed = self.run_shell(TEST_MODE="fail")
        self.assertNotEqual(failed.returncode, 0)
        self.assertFalse(self.dosage.exists())
        self.assertEqual(self.erepo.read_text(), EREPO)
        self.assertNotIn("UPDATE gene_dosage_sensitivity", failed.stdout)
        retried = self.run_shell()
        self.assertEqual(retried.returncode, 0, retried.stderr)
        self.assertIn("https://selected.invalid/GRCh38.tsv", retried.stdout)
        self.assertEqual(self.dosage.read_text(), DOSAGE.split("\n", 1)[1])

    def test_failed_or_invalid_refresh_preserves_prior_bytes_and_retry_replaces(self):
        for kind, target, payload in (("clingen_erepo", self.erepo, EREPO), ("clingen_dosage", self.dosage, DOSAGE)):
            for mode in ("fail", "empty", "invalid"):
                with self.subTest(kind=kind, mode=mode):
                    original = ("previous " + kind + " bytes\n").encode()
                    target.write_bytes(original)
                    command = 'priva_clingen_download --url "https://selected.invalid/' + ("classifications" if kind == "clingen_erepo" else "dosage") + '" --output "' + str(target) + '" --kind ' + kind
                    failed = self.run_shell(command, TEST_MODE=mode)
                    self.assertNotEqual(failed.returncode, 0)
                    self.assertEqual(target.read_bytes(), original)
                    retried = self.run_shell(command)
                    self.assertEqual(retried.returncode, 0, retried.stderr)
                    self.assertEqual(target.read_text(), payload)

    def test_current_tables_reuse_bytes_without_download(self):
        self.erepo.write_text(EREPO)
        self.dosage.write_text(DOSAGE.split("\n", 1)[1])
        INPUTS.source_receipt(self.dosage, "https://selected.invalid/GRCh38.tsv", record=True)
        before = self.erepo.read_bytes(), self.dosage.read_bytes()
        result = self.run_shell(TEST_MODE="fail")
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertNotIn("FETCH ", result.stdout)
        self.assertEqual((self.erepo.read_bytes(), self.dosage.read_bytes()), before)

    def test_changed_dosage_url_failure_preserves_prior_receipt_and_retry_rebinds(self):
        self.erepo.write_text(EREPO)
        self.dosage.write_text(DOSAGE.split("\n", 1)[1])
        INPUTS.source_receipt(self.dosage, "https://prior.invalid/GRCh38.tsv", record=True)
        sidecar = Path(str(self.dosage) + ".source.json")
        before = self.dosage.read_bytes(), sidecar.read_bytes()
        failed = self.run_shell(TEST_MODE="fail")
        self.assertNotEqual(failed.returncode, 0)
        self.assertEqual((self.dosage.read_bytes(), sidecar.read_bytes()), before)
        retried = self.run_shell()
        self.assertEqual(retried.returncode, 0, retried.stderr)
        self.assertTrue(INPUTS.source_receipt(self.dosage, "https://selected.invalid/GRCh38.tsv"))
        reused = self.run_shell(TEST_MODE="fail")
        self.assertEqual(reused.returncode, 0, reused.stderr)
        self.assertNotIn("FETCH ", reused.stdout)

    def test_first_source_receipt_rename_failure_is_recovered_before_manual_guard(self):
        for kind, target, other, payload, url in (
            ("erepo", self.erepo, self.dosage, DOSAGE, "https://erepo.clinicalgenome.org/evrepo/api/summary/classifications/download"),
            ("dosage", self.dosage, self.erepo, EREPO, "https://selected.invalid/GRCh38.tsv"),
        ):
            with self.subTest(kind=kind):
                for path in (target, Path(str(target) + ".source.json"), Path(str(target) + ".source.pending.json")):
                    if path.exists():
                        path.unlink()
                other.write_text(payload)
                if other == self.dosage:
                    INPUTS.source_receipt(other, "https://selected.invalid/GRCh38.tsv", record=True)
                failed = self.run_shell(r'''
mv() { [[ "${@: -1}" != "$TEST_FAILED_META" ]] || return 73; command mv "$@"; }
ClinGen_deploy "$TEST_CONFIG"
''', TEST_FAILED_META=str(target) + ".source.json")
                self.assertNotEqual(failed.returncode, 0)
                self.assertTrue(target.is_file())
                self.assertFalse(Path(str(target) + ".source.json").exists())
                self.assertTrue(Path(str(target) + ".source.pending.json").is_file())
                retried = self.run_shell()
                self.assertEqual(retried.returncode, 0, retried.stderr)
                self.assertNotIn("FETCH ", retried.stdout)
                self.assertFalse(Path(str(target) + ".source.pending.json").exists())
                self.assertTrue(INPUTS.source_receipt(target, url))

    def test_invalid_or_linked_pending_metadata_cannot_claim_manual_bytes(self):
        self.erepo.write_text(EREPO)
        self.dosage.write_text(DOSAGE)
        pending = Path(str(self.dosage) + ".source.pending.json")
        before = self.dosage.read_bytes()
        external = self.folder / "external.json"
        external.write_text("{}\n")
        for linked in (False, True):
            if linked:
                pending.unlink()
                pending.symlink_to(external)
            else:
                pending.write_text("{}\n")
            result = self.run_shell()
            self.assertNotEqual(result.returncode, 0)
            self.assertNotIn("FETCH ", result.stdout)
            self.assertEqual(self.dosage.read_bytes(), before)

    def test_unreceipted_manual_dosage_is_not_relabelled_for_custom_url(self):
        self.erepo.write_text(EREPO)
        self.dosage.write_text(DOSAGE)
        before = self.dosage.read_bytes()
        result = self.run_shell()
        self.assertNotEqual(result.returncode, 0)
        self.assertIn("no receipt for the nondefault URL", result.stderr)
        self.assertNotIn("FETCH ", result.stdout)
        self.assertEqual(self.dosage.read_bytes(), before)
        self.assertFalse(Path(str(self.dosage) + ".source.json").exists())

    def test_unreceipted_official_default_retains_legacy_selected_input_reuse(self):
        self.erepo.write_text(EREPO)
        self.dosage.write_text(DOSAGE.split("\n", 1)[1])
        config = json.loads(self.config.read_text())
        config["clingen_dosage_grch38_url"] = "https://ftp.clinicalgenome.org/ClinGen_gene_curation_list_GRCh38.tsv"
        self.config.write_text(json.dumps(config))
        result = self.run_shell(TEST_MODE="fail")
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertNotIn("FETCH ", result.stdout)
        self.assertFalse(Path(str(self.dosage) + ".source.json").exists())

    def test_grch37_never_uses_selected_grch38_url(self):
        config = json.loads(self.config.read_text())
        config["assembly"] = "GRCh37"
        self.config.write_text(json.dumps(config))
        self.erepo.write_text(EREPO)
        result = self.run_shell()
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertIn("ClinGen_gene_curation_list_GRCh37.tsv", result.stdout)
        self.assertNotIn("https://selected.invalid/GRCh38.tsv", result.stdout)

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

    def test_shell_forwards_selected_urls_only_to_existing_source_roles(self):
        config = json.loads(self.config.read_text())
        config.update(gene_mechanism_cache_dir=str(self.cache), gene_mechanism_raw_dir=str(self.cache / "raw"),
            ddg2p_mechanism_evidence=str(self.cache / "evidence.tsv"),
            gene_nonlof_mechanism_json=str(self.cache / "nonlof.json.gz"),
            gene_nonlof_mechanism_schema=str(self.folder / "schema.json"),
            hgnc_table=str(self.folder / "hgnc.tsv"),
            g2p_ddg2p_url="https://selected.invalid/g2p.csv",
            orphadata_gene_disease_url="https://selected.invalid/product6.xml",
            panelapp_panels_url="https://selected.invalid/panels/")
        self.config.write_text(json.dumps(config))
        (self.folder / "schema.json").write_text("{}\n")
        (self.folder / "hgnc.tsv").write_text("hgnc_id\tsymbol\n")
        for function, expected, absent in (
            ("gene_pathogenic_mechanism_cache_install", ("g2p-ddg2p", "orphadata-gene-disease"), ("clingen-dosage-grch38", "panelapp-panels")),
            ("gene_nonlof_mechanism_cache_install", ("g2p-ddg2p", "clingen-dosage-grch38", "panelapp-panels"), ("orphadata-gene-disease",)),
        ):
            with self.subTest(function=function):
                result = self.run_shell('''
priva_bundled_cache_reuse() { return 1; }
python() { printf 'BUILD %s\\n' "$*"; return 79; }
''' + function + ' "$TEST_CONFIG" sources-only')
                self.assertNotEqual(result.returncode, 0)
                self.assertIn("BUILD ", result.stdout)
                for key in expected:
                    self.assertIn("--" + key + "-url https://selected.invalid/", result.stdout)
                for key in absent:
                    self.assertNotIn("--" + key + "-url", result.stdout)


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
