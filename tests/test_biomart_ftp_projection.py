"""Bounded extraction-contract fixtures; all generated files use a data root."""
import gzip
import importlib.util
import io
import os
from pathlib import Path
import sqlite3
import tempfile
import unittest
from unittest.mock import patch
import xml.etree.ElementTree as ET

HERE = Path(__file__).resolve().parent
SPEC = importlib.util.spec_from_file_location("biomart_ftp_candidate", HERE.parent / "scripts/biomart_ftp_nm_to_enst.py")
C = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(C)


class Candidate(unittest.TestCase):
    def setUp(self):
        root = Path(os.environ.get("PRIVA_TEST_TMPDIR", tempfile.gettempdir())) / "priva-biomart116-ftp-tests"
        root.mkdir(parents=True, exist_ok=True)
        self.temp = tempfile.TemporaryDirectory(dir=root)
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)

    def database(self):
        db = sqlite3.connect(self.root / "fixture.sqlite")
        self.addCleanup(db.close)
        for alias, (_, columns) in C.PROJECTIONS.items():
            C.create_projection(db, alias, columns)
        # One protein-coding and one noncoding transcript with no RefSeq/SwissProt;
        # the last protein-coding transcript has RefSeq but no SwissProt.
        db.executemany("INSERT INTO transcripts VALUES (?,?,?,?,?,?)", (
            (1, 10, 100, "ENSG00000000010", "ENST00000000100", "protein_coding"),
            (2, 20, 200, "ENSG00000000020", "ENST00000000200", "lncRNA"),
            (3, 30, 300, "ENSG00000000030", "ENST00000000300", "protein_coding")))
        db.executemany("INSERT INTO main VALUES (?,?,?,?,?,?,?)", (
            (1, 10, 100, 1000, "ENSG00000000010", "ENST00000000100", "protein_coding"),
            (2, 20, 200, None, "ENSG00000000020", "ENST00000000200", "lncRNA"),
            (3, 30, 300, 3000, "ENSG00000000030", "ENST00000000300", "protein_coding")))
        # Genuine duplicate dimension rows must survive every combination.
        db.executemany("INSERT INTO hgnc VALUES (?,?,?)", ((1, 10, "A"), (2, 10, "A"), (3, 30, "C")))
        db.executemany("INSERT INTO swissprot VALUES (?,?,?)", ((1, 1000, "P00001"), (2, 1000, "P00001")))
        db.executemany("INSERT INTO refseq VALUES (?,?,?)", ((1, 100, "NM_000001"), (2, 100, "NM_000002"), (3, 300, "NM_000003")))
        db.commit()
        return db

    def config(self):
        row = ["95", C.DATASET, "Human genes (GRCh38.p14)", "Ensembl genes", "TableSet", "1", C.ASSEMBLY, "2026-04-30"]
        root = ET.Element("DatasetConfig", dataset=C.DATASET, datasetID="95", version=C.ASSEMBLY,
                          displayName=row[2], interfaces="default")
        for name, (key, table, field, display) in C.ATTRIBUTES.items():
            ET.SubElement(root, "AttributeDescription", internalName=name, key=key,
                          tableConstraint=table, field=field, displayName=display)
        for level, key in zip(("gene", "transcript", "translation"),
                              ("gene_id_1020_key", "transcript_id_1064_key", "translation_id_1068_key")):
            ET.SubElement(root, "MainTable").text = C.DATASET + "__" + level + "__main"
            ET.SubElement(root, "Key").text = key
        return row, root

    def test_mysql_null_empty_literal_backslash_and_binary_escapes_are_distinct(self):
        self.assertIsNone(C.mysql_value(b"\\N"))
        self.assertEqual(C.mysql_value(b""), b"")
        self.assertEqual(C.mysql_value(b"\\\\N"), b"\\N")
        self.assertEqual(C.mysql_value(b"a\\tb\\nc\\rd\\0e\\Zf\\bg\\\\h"), b"a\tb\nc\rd\0e\x1af\bg\\h")
        self.assertEqual(C.mysql_value(b"\\q"), b"q")
        with self.assertRaisesRegex(ValueError, "Truncated"):
            C.mysql_value(b"broken\\")

    def test_xml_all_five_exact_attributes_headers_and_assembly_are_required(self):
        row, root = self.config()
        self.assertEqual(C.parse_xml(ET.tostring(root), row)["main_keys"][-1], "translation_id_1068_key")
        for change in ("assembly", "header", "key", "duplicate"):
            with self.subTest(change=change):
                row, root = self.config()
                if change == "assembly":
                    row[6] = "GRCh37.p13"
                elif change == "header":
                    root.find("AttributeDescription").set("displayName", "different")
                elif change == "key":
                    root.findall("AttributeDescription")[3].set("key", "transcript_id_1064_key")
                else:
                    node = root.find("AttributeDescription")
                    ET.SubElement(root, "AttributeDescription", **node.attrib)
                with self.assertRaises(ValueError):
                    C.parse_xml(ET.tostring(root), row)

    def test_left_joins_retain_untranslated_transcripts_blank_fields_and_multiplicity(self):
        db = self.database()
        audit = C.relationship_audit(db)
        self.assertEqual(audit["transcripts_missing_from_translation_main"], 0)
        self.assertEqual(audit["noncoding_no_translation_transcripts"], 1)
        self.assertEqual(audit["untranslated_output_rows"], 1)
        output = self.root / "candidate.tsv"
        result = C.export_candidate(db, output, audit)
        rows = output.read_text().splitlines()[1:]
        self.assertEqual(result["rows"], 10)
        self.assertEqual(rows.count("ENSG00000000010\tENST00000000100\tA\tP00001\tNM_000001"), 4)
        self.assertEqual(rows.count("ENSG00000000010\tENST00000000100\tA\tP00001\tNM_000002"), 4)
        self.assertEqual(rows[-2], "ENSG00000000020\tENST00000000200\t\t\t")
        self.assertEqual(rows[-1], "ENSG00000000030\tENST00000000300\tC\t\tNM_000003")
        self.assertEqual(result["blank_field_counts"][C.HEADERS[3]], 2)

    def test_new_order_is_deterministic_without_deduplication(self):
        db = self.database()
        audit = C.relationship_audit(db)
        paths = (self.root / "one.tsv", self.root / "two.tsv")
        first = C.export_candidate(db, paths[0], audit)
        second = C.export_candidate(db, paths[1], audit)
        self.assertEqual(paths[0].read_bytes(), paths[1].read_bytes())
        self.assertEqual(first["sha256"], second["sha256"])
        self.assertEqual(len(paths[0].read_text().splitlines()), 11)

    def test_missing_untranslated_transcript_fails_closed(self):
        db = self.database()
        db.execute("DELETE FROM main WHERE transcript_id_1064_key=200")
        with self.assertRaisesRegex(ValueError, "omits transcripts"):
            C.relationship_audit(db)

    def test_inconsistent_gene_or_transcript_identity_fails_closed(self):
        db = self.database()
        db.execute("UPDATE main SET stable_id_1066='ENST99999999999' WHERE transcript_id_1064_key=200")
        with self.assertRaisesRegex(ValueError, "inconsistent"):
            C.relationship_audit(db)

    def test_duplicate_transcript_key_is_not_silently_collapsed(self):
        db = self.database()
        db.execute("INSERT INTO transcripts SELECT 4,gene_id_1020_key,transcript_id_1064_key,stable_id_1023,stable_id_1066,biotype_1064 FROM transcripts WHERE __row=1")
        with self.assertRaisesRegex(ValueError, "not unique"):
            C.relationship_audit(db)

    def test_orphan_dimension_key_fails_before_export(self):
        db = self.database()
        db.execute("INSERT INTO swissprot VALUES (3,99999,'P99999')")
        with self.assertRaisesRegex(ValueError, "Orphan.*swissprot"):
            C.relationship_audit(db)

    def test_untranslated_zero_key_has_no_fake_swissprot(self):
        db = self.database()
        db.execute("UPDATE main SET translation_id_1068_key=0 WHERE transcript_id_1064_key=200")
        C.relationship_audit(db)
        db.execute("INSERT INTO swissprot VALUES (3,0,'P99999')")
        with self.assertRaisesRegex(ValueError, "Untranslated.*SwissProt"):
            C.relationship_audit(db)

    def test_column_order_nulls_and_escaped_selected_fields_are_checked(self):
        raw = self.root
        (raw / "source.txt.gz").write_bytes(gzip.compress(b"10\t\\N\t\n20\tvalue\ttext\n"))
        db = sqlite3.connect(raw / "projection.sqlite")
        self.addCleanup(db.close)
        schema = [("gene_id_1020_key", "int(10)"), ("label", "varchar(40)"), ("other", "text")]
        stats = C.load_projection(db, raw, "projection", "source", ("label", "gene_id_1020_key", "other"), schema)
        self.assertEqual(db.execute("SELECT label,gene_id_1020_key,other FROM projection ORDER BY __row").fetchall(), [(None, 10, ""), ("value", 20, "text")])
        self.assertEqual(stats["sql_null_counts"], {"label": 1})
        self.assertEqual(stats["empty_string_counts"], {"other": 1})
        (raw / "malformed.txt.gz").write_bytes(gzip.compress(b"10\tonly_two_columns\n"))
        with self.assertRaisesRegex(ValueError, "column count"):
            C.load_projection(db, raw, "bad", "malformed", ("label",), schema)

    def test_corrupt_gzip_and_wrong_existing_raw_are_rejected_without_overwrite(self):
        raw = self.root / "x.gz"
        payload = gzip.compress(b"fixture")
        raw.write_bytes(payload[:-4])
        with self.assertRaises((EOFError, OSError)):
            C.verify_raw(raw, (raw.stat().st_size, None, None))
        before = raw.read_bytes()
        with patch.dict(C.PINS, {"x.gz": (100, None, None)}):
            with self.assertRaisesRegex(ValueError, "size mismatch"):
                C.acquire(self.root, "x.gz")
        self.assertEqual(raw.read_bytes(), before)

    def test_resume_requires_exact_http_range_and_retains_partial_on_rejection(self):
        partial = self.root / "example.gz.part"
        partial.write_bytes(b"keep")
        response = io.BytesIO(b"wrong-full-response")
        response.status = 200
        response.url = C.BASE_URL + "example.gz"
        response.headers = {"Content-Length": "20"}
        with patch.dict(C.PINS, {"example.gz": (20, None, None)}), patch.object(C, "urlopen", return_value=response):
            with self.assertRaisesRegex(ValueError, "resume range"):
                C.acquire(self.root, "example.gz")
        self.assertEqual(partial.read_bytes(), b"keep")

    def test_candidate_output_is_never_overwritten(self):
        db = self.database()
        audit = C.relationship_audit(db)
        path = self.root / "existing.tsv"
        path.write_text("preserve")
        with self.assertRaisesRegex(ValueError, "Preserve"):
            C.export_candidate(db, path, audit)
        self.assertEqual(path.read_text(), "preserve")


if __name__ == "__main__":
    unittest.main(verbosity=2)
