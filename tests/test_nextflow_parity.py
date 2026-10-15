"""Acceptance-tool tests: scientific mismatches must never be normalized away."""
import csv
from importlib.machinery import SourceFileLoader
import importlib.util
from pathlib import Path
import tempfile
import unittest

ROOT = Path(__file__).resolve().parents[1]
LOADER = SourceFileLoader("priva_compare", str(ROOT / "bin" / "priva-compare"))
SPEC = importlib.util.spec_from_loader(LOADER.name, LOADER)
COMPARE = importlib.util.module_from_spec(SPEC)
LOADER.exec_module(COMPARE)


class ParityTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.paths = [Path(self.temp.name) / name for name in ("ref.tsv", "new.tsv", "ref.acmg.tsv", "new.acmg.tsv")]
        self.header = COMPARE.FINAL_REQUIRED + ["Feature", "control_common_index"]
        self.final_rows = [
            ["1:20:A-G", "1", "20", "A", "G", "0.95", "1", "PM2_Supporting", "0.95", "Likely Pathogenic", "ENST1", "1.0"],
            ["1:20:A-G", "1", "20", "A", "G", "0.8", "1", "PM2_Supporting", "0.95", "Likely Pathogenic", "ENST2", "1.0"],
        ]
        self.matrix_header = [*COMPARE.CRITERIA, "variant_id", "frequency_conflict_PM2_BS1"]
        strengths = ["1" if criterion == "PM2" else "0" for criterion in COMPARE.CRITERIA]
        self.matrix_rows = [strengths + ["1:20:A-G", "False"], strengths + ["1:20:A-G", "False"]]
        for path in self.paths[:2]:
            self.write(path, self.header, self.final_rows)
        for path in self.paths[2:]:
            self.write(path, self.matrix_header, self.matrix_rows)

    @staticmethod
    def write(path, header, rows):
        with path.open("w", newline="") as handle:
            writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
            writer.writerow(header)
            writer.writerows(rows)

    def compare(self, **kwargs):
        return COMPARE.compare_outputs(*self.paths, **kwargs)

    def test_identical_duplicate_variant_transcripts_pass(self):
        self.assertTrue(self.compare()["identical"])

    def test_one_last_digit_of_sort_index_fails(self):
        self.final_rows[0][5] = "0.9500000000000001"
        self.write(self.paths[1], self.header, self.final_rows)
        self.assertFalse(self.compare()["identical"])

    def test_equivalent_numeric_spelling_is_not_silently_normalized(self):
        self.final_rows[0][5] = "0.950"
        self.write(self.paths[1], self.header, self.final_rows)
        self.assertFalse(self.compare()["identical"])

    def test_row_order_and_transcript_tie_order_fail(self):
        self.write(self.paths[1], self.header, list(reversed(self.final_rows)))
        self.assertFalse(self.compare()["identical"])

    def test_duplicate_row_loss_fails(self):
        self.write(self.paths[3], self.matrix_header, self.matrix_rows[:1])
        self.assertFalse(self.compare()["identical"])

    def test_duplicate_matrix_second_occurrence_strength_fails(self):
        self.matrix_rows[1][self.matrix_header.index("PM2")] = "2"
        self.write(self.paths[3], self.matrix_header, self.matrix_rows)
        result = self.compare()
        self.assertFalse(result["identical"])
        self.assertEqual(result["acmg"]["differences"][0]["key"], ("1:20:A-G", 2))

    def test_blank_is_not_nan(self):
        self.final_rows[0][7] = ""
        self.write(self.paths[0], self.header, self.final_rows)
        self.final_rows[0][7] = "NaN"
        self.write(self.paths[1], self.header, self.final_rows)
        self.assertFalse(self.compare()["identical"])

    def test_criteria_summary_order_is_exact(self):
        self.final_rows[0][7] = "PM2_Supporting;PP3"
        self.write(self.paths[0], self.header, self.final_rows)
        self.final_rows[0][7] = "PP3;PM2_Supporting"
        self.write(self.paths[1], self.header, self.final_rows)
        self.assertFalse(self.compare()["identical"])

    def test_audit_columns_are_checked(self):
        self.matrix_rows[0][-1] = "True"
        self.write(self.paths[3], self.matrix_header, self.matrix_rows)
        self.assertFalse(self.compare()["identical"])

    def test_duplicate_header_is_rejected(self):
        self.write(self.paths[1], self.header + ["Feature"], [self.final_rows[0] + ["ENST1"]])
        with self.assertRaises(ValueError):
            self.compare()

    def test_missing_criterion_column_is_rejected(self):
        self.write(self.paths[3], self.matrix_header[1:], [row[1:] for row in self.matrix_rows])
        with self.assertRaises(ValueError):
            self.compare()

    def test_zero_byte_requires_explicit_empty_case(self):
        for path in self.paths:
            path.write_bytes(b"")
        with self.assertRaises(ValueError):
            self.compare()
        self.assertTrue(self.compare(allow_empty=True)["identical"])

    def test_mixed_zero_byte_and_data_are_rejected(self):
        self.paths[3].write_bytes(b"")
        with self.assertRaises(ValueError):
            self.compare(allow_empty=True)

    def test_newline_encoding_does_not_change_tsv_values(self):
        self.paths[1].write_bytes(self.paths[1].read_bytes().replace(b"\n", b"\r\n"))
        result = self.compare()
        self.assertTrue(result["identical"])
        self.assertFalse(result["final"]["byte_identical"])


if __name__ == "__main__":
    unittest.main()
