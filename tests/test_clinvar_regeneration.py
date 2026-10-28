from __future__ import annotations

import math
from pathlib import Path
import sys
import tempfile
import unittest
import re
import subprocess
import pickle

import pandas as pd
import pysam


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

import clinvar_stat_variants as stats
from clinvar_annotation_contract import (
    AA_CSQ_FIELDS, GENE_CSQ_FIELDS, AF_CSQ_FIELDS, TABLE_INFO_FIELDS, validate_header, validate_cache,
)
from combine_annotations import arg_generator


def row(variant_id="v1", **updates):
    result = {
        "variant_id": variant_id, "SYMBOL": "TEST", "Gene": "ENSG00000000001",
        "Feature": "ENST00000000001", "BIOTYPE": "protein_coding",
        "CLNREVSTAT": "criteria_provided,_multiple_submitters,_no_conflicts",
        "DOMAINS": "Pfam:PF00001", "CLNSIG": "Pathogenic",
        "SpliceVault_SpliceAI_delta": float("nan"),
        "SpliceVault_top_events": float("nan"), "LoF": float("nan"),
        "NMD": float("nan"), "SpliceAI_pred_DS_DL": float("nan"),
        "Consequence": "missense_variant", "ref": "A", "alt": "G",
    }
    result.update(updates)
    return result


def check_explicit_summary(tmp_path):
    summary = tmp_path / "summary.tsv"
    summary.write_text(
        "#AlleleID\tAssembly\tName\n"
        "1\tGRCh37\tNM_000001.1(TEST):c.1A>G\n"
        "1\tGRCh38\tNM_000002.2(TEST):c.2A>G\n"
        "2\tGRCh38\tNM_000002.3(TEST):c.2A>G\n"
        "3\tGRCh38\tNM_000001.1(TEST):c.1A>G\n"
        "3\tGRCh38\tNM_000002.2(TEST):c.2A>G\n"
    )
    mane = tmp_path / "mane.tsv"
    mane.write_text(
        "RefSeq_nuc\tEnsembl_nuc\n"
        "NM_000001.1\tENST00000000001.1\n"
        "NM_000002.2\tENST00000000002.9\n"
    )
    mapping = stats.offer_clinvar_alleleid_to_transcript_map(str(summary), str(mane))
    assert mapping[1] == mapping[3] == "ENST00000000002"
    assert pd.isna(mapping[2])


def check_missing_splice_scores(score):
    result = stats.per_gene_stat(pd.DataFrame([row(
        SpliceVault_SpliceAI_delta=score, SpliceAI_pred_DS_DL=score,
        SpliceVault_top_events=pd.NA, LoF=pd.NA, DOMAINS=pd.NA,
    )]))
    assert result["putative_nmd_variants"] == set()
    assert result["large_aachange"] == set()
    assert result["small_aachange"] == {"v1"}
    assert result["domains"] == {}
    assert math.isnan(stats._number_or_nan(score))


def check_donor_loss_threshold(score, expected):
    result = stats.per_gene_stat(pd.DataFrame([row(SpliceAI_pred_DS_DL=score)]))
    assert result["putative_nmd_variants"] == expected


def check_original_categories():
    rows = [
        row("nmd", LoF="HC", Consequence="missense_variant&inframe_deletion", ref="A" * 16),
        row("escape", LoF="OS", NMD="NMD_escape_variant"),
        row("large", Consequence="inframe_deletion", ref="A" * 16),
        row("small", Consequence="inframe_deletion", ref="A" * 10),
        row("gap", Consequence="inframe_deletion", ref="A" * 13),
        row("synonymous", Consequence="synonymous_variant"),
        row("splice_large", SpliceVault_SpliceAI_delta=0.5,
            SpliceVault_top_events="Top1:ES:2:20%:Frameshift"),
        row("splice_small", Consequence="synonymous_variant", SpliceVault_SpliceAI_delta=0.3,
            SpliceVault_top_events="Top1:CD:-3:20%:In_frame"),
        row("one_star", CLNREVSTAT="criteria_provided,_single_submitter"),
        row("noncoding", BIOTYPE="lncRNA"),
    ]
    result = stats.per_gene_stat(pd.DataFrame(rows))
    assert result["putative_nmd_variants"] == {"nmd"}
    assert result["large_aachange"] == {"large", "splice_large"}
    assert result["small_aachange"] == {"escape", "small", "splice_small"}
    assert result["no_aachange"] == {"synonymous"}
    assert "gap" in result["clinvar_variants"]
    assert not {"one_star", "noncoding"} & result["clinvar_variants"]
    for category in ("putative_nmd_variants", "large_aachange", "small_aachange", "no_aachange"):
        assert result["transcripts"]["ENST00000000001"][category] == result[category]
        assert result["domains"]["ENSG00000000001:Pfam:PF00001"][category] == result[category]


def header_vcf(path, fields, *, clncsq=False, population=False):
    header = pysam.VariantHeader()
    header.contigs.add("chr1")
    for name in ("ALLELEID", "CLNSIG", "CLNREVSTAT"):
        header.add_meta("INFO", items=[("ID", name), ("Number", "."), ("Type", "String"), ("Description", name)])
    header.add_meta("INFO", items=[("ID", "CSQ"), ("Number", "."), ("Type", "String"),
                                  ("Description", "VEP annotations. Format: " + "|".join(sorted(fields)))])
    if clncsq:
        header.add_meta("INFO", items=[("ID", "CLNCSQ"), ("Number", "."), ("Type", "String"),
                                      ("Description", "ClinVar annotations. Format: Feature|SYMBOL")])
    if population:
        for name in TABLE_INFO_FIELDS:
            header.add_meta("INFO", items=[("ID", name), ("Number", "."), ("Type", "Float"), ("Description", name)])
    with pysam.VariantFile(str(path), "w", header=header):
        pass


def check_clinvar_csq(test, tmp_path):
    path = tmp_path / "source.vcf"
    header_vcf(path, {"Feature", "SYMBOL"})
    with test.assertRaises(KeyError):
        list(arg_generator(path, 1, {}, {}))
    assert list(arg_generator(path, 1, {}, {}, require_clinvar_csq=False)) == []
    header_vcf(path, {"Feature", "SYMBOL"}, clncsq=True)
    assert list(arg_generator(path, 1, {}, {})) == []


def check_annotation_contract(test, tmp_path):
    path = tmp_path / "source.vcf"
    header_vcf(path, AA_CSQ_FIELDS | {"SYMBOL"})
    assert validate_header(path, "aa")["success"]
    with test.assertRaisesRegex(ValueError, "SpliceVault_SpliceAI_delta"):
        validate_header(path, "vep")
    header_vcf(path, AA_CSQ_FIELDS | GENE_CSQ_FIELDS)
    assert validate_header(path, "vep")["success"]
    with test.assertRaisesRegex(ValueError, "AF_joint"):
        validate_header(path, "complete")
    header_vcf(path, AA_CSQ_FIELDS | GENE_CSQ_FIELDS | AF_CSQ_FIELDS, population=True)
    assert validate_header(path, "complete")["success"]


class ClinVarRegenerationTests(unittest.TestCase):
    def setUp(self):
        was_disabled = stats.logger.disabled
        stats.logger.disabled = True
        self.addCleanup(setattr, stats.logger, "disabled", was_disabled)
        self.temporary = tempfile.TemporaryDirectory()
        self.addCleanup(self.temporary.cleanup)
        self.directory = Path(self.temporary.name)

    def test_explicit_summary_keeps_grch38_exact_version_mapping_and_last_row(self):
        check_explicit_summary(self.directory)

    def test_missing_splice_scores_remain_absent(self):
        for score in (pd.NA, None, float("nan")):
            with self.subTest(score=str(score)):
                check_missing_splice_scores(score)

    def test_donor_loss_threshold_keeps_strict_greater_than(self):
        for score, expected in ((0.5, set()), (0.50001, {"v1"})):
            with self.subTest(score=score):
                check_donor_loss_threshold(score, expected)

    def test_original_category_precedence_and_review_filtering(self):
        check_original_categories()

    def test_original_nmd_escaping_label_behavior_is_preserved(self):
        result = stats.per_gene_stat(pd.DataFrame([row(LoF="HC", NMD="NMD_escaping_variant")]))
        self.assertEqual(result["putative_nmd_variants"], {"v1"})

    def test_full_annotation_uses_one_native_glob_directory_list_and_same_plugin_order(self):
        source = (ROOT / "scripts/annotation_vcf.sh").read_text()
        function = source.split("function anno_VEP_data() {", 1)[1].split("function Calculate_CADD", 1)[0]
        command = function.split('vep -i "${input_vcf}"', 1)[1]
        line = next(value.strip().rstrip("\\").strip() for value in command.splitlines()
                    if value.strip().startswith("--dir_plugins "))
        shell = ('vep_plugins_dir=/configured/plugins; loftee_repo=/configured/loftee; '
                 'capture() { printf "%s\\0" "$@"; }; capture ' + line)
        result = subprocess.run(["bash", "-c", shell], check=True, capture_output=True)
        self.assertEqual(result.stdout.split(b"\0")[:-1],
                         [b"--dir_plugins", b"/configured/plugins /configured/loftee"])
        glob = subprocess.run(["perl", "-e", 'print join("\\0",glob(shift @ARGV));',
                               "/configured/plugins /configured/loftee"],
                              check=True, capture_output=True)
        self.assertEqual(glob.stdout.split(b"\0"), [b"/configured/plugins", b"/configured/loftee"])
        self.assertEqual(re.findall(r"^\s+-plugin (\w+)", command, re.MULTILINE),
                         ["UTRAnnotator", "LOEUF", "AlphaMissense", "LoF", "SpliceAI",
                          "PrimateAI", "SpliceVault", "Conservation", "NMD"])

    def test_clinvar_csq_default_stays_strict_and_source_mode_is_explicit(self):
        check_clinvar_csq(self, self.directory)

    def test_annotation_contract_rejects_basic_vep_and_distinguishes_af_stage(self):
        check_annotation_contract(self, self.directory)

    def test_cache_validation_rejects_partial_empty_and_wrong_stage_outputs(self):
        path = self.directory / 'stats.pkl'
        path.write_bytes(b'\x80\x04')
        with self.assertRaises(EOFError):
            validate_cache(path, 'aa')
        path.write_bytes(pickle.dumps({}))
        with self.assertRaises(ValueError):
            validate_cache(path, 'aa')
        path.write_bytes(pickle.dumps({'ENST1': {'1': {'p.A1V': {'CLNSIG': ['Pathogenic'], 'CLNREVSTAT': ['reviewed_by_expert_panel']}}}}))
        self.assertTrue(validate_cache(path, 'aa')['success'])
        with self.assertRaises(ValueError):
            validate_cache(path, 'splice')
        path.write_bytes(pickle.dumps({'ENSG1': None}))
        self.assertTrue(validate_cache(path, 'gene-af')['success'])


if __name__ == "__main__":
    unittest.main()
