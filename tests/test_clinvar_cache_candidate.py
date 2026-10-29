"""Bounded ClinVar installation candidates stay separate from runtime annotation."""
from __future__ import annotations

import json
import os
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest

import pysam
import yaml

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
import clinvar_gene_refresh as refresh
import clinvar_annotation_contract as contract
import spliceai_cache_input as spliceai


def write_vcf(path, *, annotated=False, version2=False):
    path.parent.mkdir(parents=True, exist_ok=True)
    text = "##fileformat=VCFv4.2\n##fileDate=20260430\n##contig=<ID=chr1,length=1000>\n"
    if version2:
        text += "##SpliceAI2Version=2\n##assembly=GRCh38\n"
    if annotated:
        fields = sorted(contract.AA_CSQ_FIELDS | contract.GENE_CSQ_FIELDS)
        text += '##INFO=<ID=CSQ,Number=.,Type=String,Description="Format: ' + "|".join(fields) + '">\n'
        for name in ("CLNSIG", "CLNREVSTAT", "ALLELEID"):
            text += f'##INFO=<ID={name},Number=1,Type=String,Description="fixture">\n'
    text += "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\nchr1\t10\t.\tA\tC\t.\t.\t.\n"
    with pysam.BGZFile(str(path), "wb") as output:
        output.write(text.encode())
    pysam.tabix_index(str(path), preset="vcf", force=True)


class ClinVarCacheCandidateTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory(prefix="clinvar_candidate_",
            dir=os.environ.get("PRIVA_TEST_TMPDIR", os.environ.get("TMPDIR")))
        self.addCleanup(self.temp.cleanup)
        self.folder = Path(self.temp.name)
        self.cache = self.folder / "cache"
        self.vcf = self.cache / "clinvar/clinvar.hg38.vep.vcf.gz"
        write_vcf(self.vcf)
        self.snv = self.cache / "spliceai/snv.vcf.gz"
        self.indel = self.cache / "spliceai/indel.vcf.gz"
        write_vcf(self.snv)
        write_vcf(self.indel)
        self.config = self.folder / "config.yaml"
        self.settings = {
            "assembly": "hg38", "vep_installed_version": "116", "cache_root": str(self.cache),
            "clinvar_vcf_dir": str(self.vcf.parent), "clinvar_vcf": str(self.vcf), "threads": 1,
            "spliceai_snv_prescore": str(self.snv), "spliceai_indel_prescore": str(self.indel),
            "ref_genome": str(self.cache / "reference.fa"), "vep_cache_dir": str(self.cache / "vep"),
            "vep_plugins_dir": str(self.cache / "plugins"), "vep_plugins_cachedir": str(self.cache / "plugin-cache"),
            "loftee_repo": str(self.cache / "loftee"), "loeuf_prescore": "loeuf.tsv.gz",
            "pext_enabled": False,
        }
        self.save_config()
        self.trace = self.folder / "vep-argv.json"
        self.annotated = self.folder / "vep-fixture.vcf.gz"
        write_vcf(self.annotated, annotated=True)

    def save_config(self):
        self.config.write_text(yaml.safe_dump(self.settings, sort_keys=False))

    def v2(self):
        outputs = {}
        for role, path in (("snv", self.snv), ("indel", self.indel)):
            write_vcf(path, version2=True)
            outputs[role] = {"path": str(path), "bytes": path.stat().st_size,
                "sha256": refresh.digest(path), "index_path": str(path) + ".tbi",
                "index_sha256": refresh.digest(Path(str(path) + ".tbi"))}
        receipt = {"schema_version": 1, "resource": "SpliceAI2", "version": "2.0", "assembly": "GRCh38",
            "status": "prepared_annotation_only", "coverage": "complete_published", "scientific_qualified": False,
            "source_revision": spliceai.REVISION, "source_manifest_sha256": "a" * 64, "outputs": outputs}
        (self.snv.parent / "adapter.ready.json").write_text(json.dumps(receipt))

    def begin(self):
        return refresh.candidate_begin(self.config, ROOT, self.vcf)

    def shell(self, function="anno_VEP_cache_candidate", extra=(), env=None):
        # Scientific executable is a tiny deterministic fixture. Sort/tabix run
        # for real; all original VEP argv are retained for comparison.
        body = r'''
load_script() { source "$1"; }; load_script "$1" >/dev/null || exit
check_path() { return 0; }
check_vcf_infotags() { [[ "${CANDIDATE_TAGS_ALREADY_PRESENT:-0}" == 1 ]]; }
randomID() { printf 'FIXED\n'; }
announce_remove_tmps() { rm -f -- "$@"; }
display_vcf() { return 0; }
vep() {
    python - "$@" <<'PY'
import gzip, json, os, sys
args = sys.argv[1:]
with open(os.environ['CANDIDATE_TRACE'], 'w') as handle: json.dump(args, handle)
with gzip.open(os.environ['CANDIDATE_ANNOTATED'], 'rb') as source, open(args[args.index('-o') + 1], 'wb') as output:
    output.write(source.read())
PY
}
"$2" --config "$CANDIDATE_CONFIG" --input_vcf "$CANDIDATE_INPUT" --threads 1 "${@:3}"
'''
        environment = dict(os.environ, CANDIDATE_CONFIG=str(self.config), CANDIDATE_INPUT=str(self.vcf),
                           CANDIDATE_TRACE=str(self.trace), CANDIDATE_ANNOTATED=str(self.annotated))
        if env:
            environment.update(env)
        return subprocess.run(["bash", "-c", body, "candidate-fixture", str(ROOT / "scripts/annotation_vcf.sh"),
                               function, *extra], env=environment, capture_output=True, text=True, timeout=30)

    def test_source_probe_requires_no_config_or_cache_paths(self):
        result = refresh.candidate_route(ROOT)
        self.assertEqual(result["status"], "preparation_route_available_inputs_not_checked")
        self.assertIs(result["scientific_qualified"], False)
        self.assertEqual(result["runtime_guard_sha256"], "7adeabcab3037d248cc1460fd065bc5220158a4c5d4e9f106065568b925b5eac")
        self.assertEqual(result["helper_source"], refresh.identity(ROOT / "scripts/clinvar_gene_refresh.py"))

    def test_profile_binds_loeuf_path_actually_selected_by_original_vep_shell(self):
        self.settings["loeuf_prescore"] = str(self.cache / "ignored-directory/loeuf.tsv.gz")
        self.save_config()
        context = refresh.candidate_context(self.config, ROOT, self.vcf)
        files = context["scientific_profile"]["files"]
        actual = str(ROOT / "data/loeuf/loeuf.tsv.gz")
        self.assertIn(actual, files)
        self.assertNotIn(self.settings["loeuf_prescore"], files)
        self.assertNotIn("sha256", files[actual])

    def test_managed_v2_candidate_finishes_and_runtime_cannot_rewrite_its_receipt(self):
        self.v2()
        config_before = self.config.read_bytes()
        result = self.shell(env={"CANDIDATE_TAGS_ALREADY_PRESENT": "1"})
        self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
        receipt = json.loads(refresh.candidate_receipt_path(self.vcf).read_text())
        self.assertEqual(receipt["status"], "prepared_annotation_only")
        self.assertEqual(receipt["stage"], "vep_before_gnomad_enrichment")
        self.assertEqual(receipt["spliceai_format"], "vcf_spliceai2")
        self.assertIs(receipt["scientific_qualified"], False)
        self.assertEqual(receipt["output"], refresh.identity(self.vcf))
        self.assertFalse(refresh.candidate_receipt_path(self.vcf, pending=True).exists())
        self.assertEqual(self.config.read_bytes(), config_before)
        # Existing CSQ/mtime shortcuts are disabled for candidates, so VEP did
        # execute even though the fixture's tag check says already annotated.
        self.assertTrue(self.trace.exists())
        self.trace.unlink()
        ready = refresh.candidate_receipt_path(self.vcf).read_bytes()
        # Runtime accepts SpliceAI2 with its own cutoffs; an exported purpose
        # still cannot turn it into a candidate preparation.
        runtime = self.shell("anno_VEP_data", env={"priva_vep_purpose": "cache_candidate"})
        self.assertEqual(runtime.returncode, 0, runtime.stdout + runtime.stderr)
        self.assertIn("SpliceAI2 cutoffs", runtime.stderr)
        self.assertEqual(refresh.candidate_receipt_path(self.vcf).read_bytes(), ready)
        self.assertFalse(refresh.candidate_receipt_path(self.vcf, pending=True).exists())

    def test_candidate_uses_original_vep_arguments_and_cutoff(self):
        result = self.shell()
        self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
        candidate_args = json.loads(self.trace.read_text())
        self.trace.unlink()
        runtime = self.shell("anno_VEP_data")
        self.assertEqual(runtime.returncode, 0, runtime.stdout + runtime.stderr)
        self.assertEqual(json.loads(self.trace.read_text()), candidate_args)
        self.assertEqual([item for item in candidate_args if item.startswith("SpliceAI,")],
                         [f"SpliceAI,snv={self.snv},indel={self.indel},cutoff=0.5"])
        plugins = [candidate_args[index + 1] for index, value in enumerate(candidate_args[:-1]) if value == "-plugin"]
        self.assertEqual([value.split(",", 1)[0] for value in plugins],
                         ["UTRAnnotator", "LOEUF", "AlphaMissense", "LoF", "SpliceAI", "PrimateAI", "SpliceVault", "Conservation", "NMD", "MaveDB"])

    def test_private_entry_cannot_take_candidate_mode_from_environment(self):
        self.v2()
        result = self.shell("_priva_anno_VEP_data", env={"priva_vep_purpose": "cache_candidate"})
        self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
        self.assertIn("SpliceAI2 cutoffs", result.stderr)
        self.assertFalse(refresh.candidate_receipt_path(self.vcf).exists())
        self.assertFalse(refresh.candidate_receipt_path(self.vcf, pending=True).exists())

    def test_candidate_rejects_scientific_cli_overrides(self):
        result = self.shell(extra=("--spliceai_snv_prescore", str(self.snv)))
        self.assertNotEqual(result.returncode, 0)
        self.assertIn("Unsupported ClinVar candidate option", result.stdout + result.stderr)
        self.assertFalse(self.trace.exists())
        self.assertFalse(refresh.candidate_receipt_path(self.vcf, pending=True).exists())

    def test_v2_needs_completed_same_pair_adapter_receipt(self):
        self.v2()
        (self.snv.parent / "adapter.ready.json").unlink()
        with self.assertRaises((OSError, ValueError)):
            self.begin()
        self.assertFalse(refresh.candidate_receipt_path(self.vcf, pending=True).exists())

    def test_candidate_cannot_overwrite_an_external_or_unselected_vcf(self):
        outside = self.folder / "outside.vcf.gz"
        write_vcf(outside)
        with self.assertRaisesRegex(ValueError, "selected managed VCF"):
            refresh.candidate_begin(self.config, ROOT, outside)
        self.settings["clinvar_vcf"] = str(outside)
        self.save_config()
        with self.assertRaisesRegex(ValueError, "selected managed VCF"):
            refresh.candidate_begin(self.config, ROOT, outside)

    def test_symlink_output_escape_is_rejected(self):
        outside = self.folder / "outside.vcf.gz"
        write_vcf(outside)
        self.vcf.unlink()
        self.vcf.symlink_to(outside)
        with self.assertRaisesRegex(ValueError, "selected managed VCF"):
            self.begin()

    def test_changed_config_or_adapter_cannot_publish_completed_candidate(self):
        for change in ("config", "adapter"):
            with self.subTest(change=change):
                self.v2()
                self.begin()
                write_vcf(self.vcf, annotated=True)
                if change == "config":
                    self.settings["threads"] = 2
                    self.save_config()
                else:
                    path = self.snv.parent / "adapter.ready.json"
                    value = json.loads(path.read_text())
                    value["fixture_changed"] = True
                    path.write_text(json.dumps(value))
                with self.assertRaisesRegex(ValueError, "changed during preparation"):
                    refresh.candidate_complete(self.config, ROOT, self.vcf)
                self.assertFalse(refresh.candidate_receipt_path(self.vcf).exists())
                self.assertTrue(refresh.candidate_receipt_path(self.vcf, pending=True).exists())

    def test_incomplete_annotation_has_no_ready_receipt(self):
        self.begin()
        with self.assertRaisesRegex(ValueError, "CSQ"):
            refresh.candidate_complete(self.config, ROOT, self.vcf)
        self.assertFalse(refresh.candidate_receipt_path(self.vcf).exists())


if __name__ == "__main__":
    unittest.main()
