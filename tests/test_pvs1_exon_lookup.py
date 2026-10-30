from __future__ import annotations

from pathlib import Path
import sys
import unittest


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

from acmg_pvs1_null_variant import (  # noqa: E402
    downstream_exon_patho_af,
    exon_patho_af_entry,
    span_functional_domains,
)


# Exons 1, 34, 35 and 36 of a 36-exon transcript, all with common pathogenic AFs
# (median 0.05, above the 0.01 threshold), keyed the way the cache stores VEP's EXON text.
CACHE = {
    "ENST1": {
        "1/36": (0.05, 0.05, 0.05, 0.05, (0.05,)),
        "34/36": (0.05, 0.05, 0.05, 0.05, (0.05,)),
        "35/36": (0.05, 0.05, 0.05, 0.05, (0.05,)),
        "36/36": (0.05, 0.05, 0.05, 0.05, (0.05,)),
        "35-36/36": (0.05, 0.05, 0.05, 0.05, (0.05,)),
    }
}


def row(**values):
    base = {
        "Consequence": "missense_variant", "Feature": "ENST1", "Gene": "ENSG1", "chrom": "chr1", "pos": 1,
        "EXON": "34/36", "Protein_position": "10", "DOMAINS": float("nan"),
        "5UTR_lof": False, "5UTR_frameshift": False, "5UTR_len_changing": False, "5UTR_span_intol_domain": False,
        "splicing_lof": False, "splicing_frameshift": False, "splicing_len_changing": False,
        "splicing_span_intol_domain": False, "splicing_affected_exons": "",
    }
    base.update(values)
    return base


class ExonLookupTest(unittest.TestCase):
    def test_an_exon_number_finds_the_cache_key_of_its_transcript(self) -> None:
        self.assertEqual(exon_patho_af_entry(CACHE, "ENST1", 34)[0], 0.05)
        self.assertEqual(exon_patho_af_entry(CACHE, "ENST1", "1")[0], 0.05)
        self.assertEqual(exon_patho_af_entry(CACHE, "ENST1", 2), (0, ))
        self.assertEqual(exon_patho_af_entry(CACHE, "ENST_absent", 1), (0, ))

    def test_downstream_exons_use_their_cache_entries(self) -> None:
        # Exons 34 to 36 all have common pathogenic AFs, so none counts as rare.
        self.assertFalse(downstream_exon_patho_af(row(EXON="34/36"), CACHE, logic="any", threshold=0.01))
        self.assertFalse(downstream_exon_patho_af(row(EXON="34-35/36"), CACHE, logic="any", threshold=0.01))
        # Exon 33 is not in the cache: its default AF 0 counts as rare, as before.
        self.assertTrue(downstream_exon_patho_af(row(EXON="33/36"), CACHE, logic="any", threshold=0.01))

    def test_five_prime_utr_variants_read_exon_one(self) -> None:
        _, exon_frequent_patho = span_functional_domains(row(**{"5UTR_lof": True}),
                                                         clinvar_patho_exon_af_dict=CACHE)
        self.assertFalse(exon_frequent_patho)

    def test_splicing_variants_read_the_skipped_exons(self) -> None:
        _, exon_frequent_patho = span_functional_domains(row(splicing_lof=True, splicing_affected_exons="35,36"),
                                                         clinvar_patho_exon_af_dict=CACHE)
        self.assertFalse(exon_frequent_patho)

    def test_lookups_by_vep_exon_text_are_unchanged(self) -> None:
        _, exon_frequent_patho = span_functional_domains(row(EXON="35-36/36"), clinvar_patho_exon_af_dict=CACHE)
        self.assertFalse(exon_frequent_patho)
        _, exon_frequent_patho = span_functional_domains(row(EXON="2/36"), clinvar_patho_exon_af_dict=CACHE)
        self.assertTrue(exon_frequent_patho)


if __name__ == "__main__":
    unittest.main()
