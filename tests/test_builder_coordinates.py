from pathlib import Path
import shutil
import sys
import tempfile
import unittest
from unittest import mock

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

import SpliceVault_tsv_vcf_conversion as splicevault  # noqa: E402
import build_gofcards_exact_gof_cache as gofcards  # noqa: E402


def identity_liftover(input_vcf, output_vcf, input_bed, output_bed, chain_file, reference_fasta):
    shutil.copyfile(input_vcf, output_vcf)
    shutil.copyfile(input_bed, output_bed)


class SpliceVaultSitesPerRow(unittest.TestCase):
    def test_each_event_keeps_its_own_lifted_site(self):
        chunk = pd.DataFrame({
            "#chrom": ["chr1", "chr1", "chr1"], "pos": [1000, 1000, 2000],
            "ref": ["A", "A", "C"], "alt": ["G", "G", "T"],
            "splicevault_site_loc": ["chr1:1100", "chr1:1200", "chr1:2100"],
        })
        with tempfile.TemporaryDirectory() as temp, \
                mock.patch.object(splicevault, "liftover_coordinates", identity_liftover):
            updated = splicevault.process_and_update_chunk(chunk.copy(), "chain", "fasta", temp)
        self.assertEqual(list(updated["splicevault_site_loc"]), ["chr1:1100", "chr1:1200", "chr1:2100"])
        self.assertEqual(list(updated.columns), list(chunk.columns))


class GoFCardsEndpointCoordinate(unittest.TestCase):
    def test_a_reference_checked_endpoint_position_keys_the_grch38_genomic_route(self):
        row = {
            "HGNC_Symbol": "GENE1", "allele_key": "SNV|1|100|A|G", "vcf_id": "v1",
            "hg19_chrom": "1", "source_hg19_pos": "100", "source_hg19_ref": "A", "source_hg19_alt": "G",
            "match_eligibility": "eligible", "gene_match_status": "match", "reject_reason": "",
            "VEP_HGNC_Symbol": "GENE1", "gofcards_accession_id": "", "VEP_assembly": "",
            "hg38_chrom": "1", "hg38_pos": "5100", "hg38_ref": "A", "hg38_alt": "G",
            "liftover_status": "unmapped", "hg38_refalt_status": "endpoint_ref_match",
        }
        cache = gofcards._nest(pd.DataFrame([row]), {}, {}, gofcards.HgncResolver(None), {}, {"v1"})
        variant = next(iter(cache["genes"]["GENE1"]["variants"].values()))
        self.assertEqual(variant["record"]["liftover_status"], "unmapped")
        self.assertEqual(variant["assemblies"]["hg38"]["genomic"],
                         {"chrom": "1", "pos": 5100, "ref": "A", "alt": "G", "status": "endpoint_ref_match"})


if __name__ == "__main__":
    unittest.main()
