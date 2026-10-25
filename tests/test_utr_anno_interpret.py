"""Paired UTR effects must retain their individual interpretation policies."""
import importlib.util
from pathlib import Path
import unittest

import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
spec = importlib.util.spec_from_file_location("utr_interpret_test", ROOT / "scripts/utr_anno_interpret.py")
utr = importlib.util.module_from_spec(spec)
spec.loader.exec_module(utr)
utr.logger.disabled = True
START = "5_prime_UTR_premature_start_codon_gain_variant"
STOP = "5_prime_UTR_uORF_stop_codon_loss_variant"
SHIFT = "5_prime_UTR_uORF_frameshift_variant"
WEAK_START = "KozakStrength=Weak:type=OutOfFrame_oORF"
STRONG_START = "KozakStrength=Strong:type=OutOfFrame_oORF"
STRONG_STOP = "KozakStrength=Strong:AltStop=False:FrameWithCDS=outOfFrame"
STRONG_SHIFT = "KozakStrength=Strong:alt_type=OutOfFrame_oORF:alt_type_length=NA"


class PairedEffectsTests(unittest.TestCase):
    def interpret(self, effects):
        row = pd.Series({"5UTR_consequence": "&".join(e[0] for e in effects),
                         "5UTR_annotation": "&".join(e[1] for e in effects)})
        return utr.interpret_utr_annotation(row)

    def test_weak_start_does_not_hide_strong_stop_loss_or_frameshift(self):
        for label, details in ((STOP, STRONG_STOP), (SHIFT, STRONG_SHIFT)):
            for effects in (((START, WEAK_START), (label, details)),
                            ((label, details), (START, WEAK_START))):
                with self.subTest(effects=effects):
                    self.assertEqual(self.interpret(effects), (True, True, True, False))

    def test_repeated_label_preserves_multiple_orfs(self):
        self.assertEqual(self.interpret(((START, WEAK_START), (START, STRONG_START))),
                         (True, True, True, False))

    def test_single_effect_rules_and_uninterpreted_effects_unchanged(self):
        for label, details, expected in (
            (START, WEAK_START, False), (START, STRONG_START, True),
            (STOP, STRONG_STOP, True), (SHIFT, STRONG_SHIFT, True),
            ("5_prime_UTR_uORF_stop_codon_gain_variant", STRONG_STOP, False),
            ("5_prime_UTR_start_codon_loss_variant", STRONG_START, False),
            ("-", "-", False),
        ):
            with self.subTest(label=label, details=details):
                self.assertEqual(self.interpret(((label, details),)),
                                 (expected, expected, expected, False))


if __name__ == "__main__":
    unittest.main()
