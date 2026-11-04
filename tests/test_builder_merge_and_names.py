from pathlib import Path
import sys
import unittest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

from inject_clinvar_into_gofcards import clinvar_block  # noqa: E402
from stat_gene_patho_afs import merge_results  # noqa: E402


class GenePathogenicAfMerge(unittest.TestCase):
    def test_none_marks_a_gene_without_a_nonzero_af(self):
        high, low = {"af": 0.3}, {"af": 0.1}
        self.assertEqual(merge_results([{"A": low}, {"A": None}]), {"A": low})
        self.assertEqual(merge_results([{"A": None}, {"A": high}]), {"A": high})
        self.assertEqual(merge_results([{"A": None}, {"A": None}]), {"A": None})
        self.assertEqual(merge_results([{"A": low}, {"A": high}]), {"A": high})
        self.assertEqual(merge_results([{"A": high}, {"A": low}]), {"A": high})


class ClinvarBlockName(unittest.TestCase):
    def test_variation_name_comes_from_the_parsed_record(self):
        payload = {"variation": {"vcv_accession": "VCV000001", "variation_id": "1",
                                 "variation_name": "NM_000001.1(GENE):c.1A>G (p.Met1Val)", "hgvs": []},
                   "condition_assertions": []}
        self.assertEqual(clinvar_block(payload)["variation_name"], "NM_000001.1(GENE):c.1A>G (p.Met1Val)")


if __name__ == "__main__":
    unittest.main()
