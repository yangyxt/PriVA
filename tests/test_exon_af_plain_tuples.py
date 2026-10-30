from __future__ import annotations

from pathlib import Path
import pickle
import pickletools
import sys
import unittest


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

from stat_exon_patho_afs import calculate_stats  # noqa: E402


class ExonAfStatisticsTest(unittest.TestCase):
    def test_statistics_keep_their_values_and_store_the_afs_as_plain_floats(self) -> None:
        stats = calculate_stats({"ENST1": {"4/36": [0.5, 0.25, 0.125], "5/36": [0.01]}, "ENST2": {"1/2": []}})

        self.assertEqual(stats["ENST1"]["4/36"], (0.25, 0.2916666666666667, 0.5, 0.125, (0.5, 0.25, 0.125)))
        self.assertEqual(stats["ENST1"]["5/36"], (0.01, 0.01, 0.01, 0.01, (0.01,)))
        self.assertEqual(stats["ENST2"], {})
        for value in stats["ENST1"].values():
            self.assertIsInstance(value[4], tuple)
            self.assertTrue(all(type(af) is float for af in value[4]))

    def test_the_pickle_refers_to_no_numpy_object(self) -> None:
        stats = calculate_stats({"ENST1": {"4/36": [0.5, 0.25]}})
        opcodes = {opcode.name for opcode, _, _ in pickletools.genops(pickle.dumps(stats, protocol=4))}

        self.assertNotIn("GLOBAL", opcodes)
        self.assertNotIn("STACK_GLOBAL", opcodes)


if __name__ == "__main__":
    unittest.main()
