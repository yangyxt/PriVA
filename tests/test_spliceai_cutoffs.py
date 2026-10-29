"""SpliceAI v1 and SpliceAI2 scores each use their own cutoffs in PriVA's rules."""
import gzip
import os
from pathlib import Path
import sys
import tempfile
import unittest

import pandas as pd
import yaml

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
from spliceai_cutoffs import CUTOFFS, cutoffs_for_config
from acmg_pp3_bp4_bp7_insilico import PP3_BP4_criteria
import clinvar_stat_variants as stats

HEADER = b"##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
V2_HEADER = HEADER.replace(b"#CHROM", b"##SpliceAI2Version=2\n#CHROM")


class SpliceAICutoffTests(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory(prefix="spliceai_cutoffs_",
            dir=os.environ.get("PRIVA_TEST_TMPDIR", os.environ.get("TMPDIR")))
        self.addCleanup(self.temporary.cleanup)
        self.folder = Path(self.temporary.name)

    def config(self, snv_header, indel_header, assembly="hg38"):
        snv, indel = self.folder / "snv.vcf.gz", self.folder / "indel.vcf.gz"
        snv.write_bytes(gzip.compress(snv_header))
        indel.write_bytes(gzip.compress(indel_header))
        path = self.folder / "config.yaml"
        path.write_text(yaml.safe_dump({"assembly": assembly, "spliceai_snv_prescore": str(snv),
                                        "spliceai_indel_prescore": str(indel)}))
        return path

    def test_each_version_selects_its_own_cutoffs(self):
        self.assertEqual(cutoffs_for_config(self.config(HEADER, HEADER)), CUTOFFS["vcf_v1"])
        self.assertEqual(cutoffs_for_config(self.config(V2_HEADER, V2_HEADER)), CUTOFFS["vcf_spliceai2"])

    def test_original_cutoffs_apply_without_score_files(self):
        self.assertEqual(cutoffs_for_config(None), CUTOFFS["vcf_v1"])
        self.assertEqual(cutoffs_for_config(self.folder / "missing.yaml"), CUTOFFS["vcf_v1"])
        self.assertEqual(CUTOFFS["vcf_v1"], {"splice_altering": 0.8, "nearby_acceptor_rescue": 0.5,
                                             "known_pathogenic_comparison": 0.5, "no_splicing_impact": 0.1,
                                             "donor_loss_putative_nmd": 0.5})

    def test_mixed_versions_and_wrong_assembly_are_refused(self):
        with self.assertRaisesRegex(ValueError, "different versions"):
            cutoffs_for_config(self.config(V2_HEADER, HEADER))
        with self.assertRaisesRegex(ValueError, "configured assembly is GRCh37"):
            cutoffs_for_config(self.config(V2_HEADER, V2_HEADER, assembly="hg19"))

    def test_spliceai2_uses_the_authors_thresholds_for_every_rule(self):
        # SpliceAI2 authors: 0.1 high recall, 0.25 balance, 0.5 high precision
        # (in place of SpliceAI v1's 0.2, 0.5 and 0.8).
        self.assertEqual(CUTOFFS["vcf_spliceai2"], {"splice_altering": 0.5, "nearby_acceptor_rescue": 0.25,
                                                    "known_pathogenic_comparison": 0.25, "no_splicing_impact": 0.1,
                                                    "donor_loss_putative_nmd": 0.25})

    def test_no_splicing_impact_uses_the_given_cutoff(self):
        # One splice-region variant whose CADD score is neither low nor high, so
        # only SpliceAI's largest delta score (0.07) decides BP4.
        frame = pd.DataFrame({
            "Consequence": ["splice_region_variant&intron_variant"], "LoF": [""],
            "PrimateAI": [0.0], "am_pathogenicity": [0.0], "CADD_phred": [18.0], "CADD_reg_phred": [10.0],
            "am_class": [""], "vep_consq_lof": [False], "splicing_lof": [False], "5UTR_lof": [False],
            "CLNSIG": [""], "CLNREVSTAT": [""], "SpliceAI_pred_DS_AG": [0.07], "SpliceAI_pred_DS_AL": [0.0],
            "SpliceAI_pred_DS_DG": [0.01], "SpliceAI_pred_DS_DL": [0.0],
        })
        self.assertEqual(PP3_BP4_criteria(frame)[1].tolist(), [1])
        self.assertEqual(PP3_BP4_criteria(frame, spliceai_no_impact_cutoff=0.07)[1].tolist(), [0])

    def test_donor_loss_counted_as_putative_nmd_uses_the_given_cutoff(self):
        row = {"variant_id": "v1", "SYMBOL": "TEST", "Gene": "ENSG00000000001", "Feature": "ENST00000000001",
               "BIOTYPE": "protein_coding", "CLNREVSTAT": "criteria_provided,_multiple_submitters,_no_conflicts",
               "DOMAINS": "Pfam:PF00001", "CLNSIG": "Pathogenic", "SpliceVault_SpliceAI_delta": float("nan"),
               "SpliceVault_top_events": float("nan"), "LoF": float("nan"), "NMD": float("nan"),
               "SpliceAI_pred_DS_DL": 0.4, "Consequence": "missense_variant", "ref": "A", "alt": "G"}
        self.assertEqual(stats.per_gene_stat(pd.DataFrame([row]))["putative_nmd_variants"], set())
        self.assertEqual(stats.per_gene_stat(pd.DataFrame([row]), 0.35)["putative_nmd_variants"], {"v1"})


if __name__ == "__main__":
    unittest.main()
