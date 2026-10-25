"""Real bounded BGZF/TBI repair and original-command PrimateAI parity fixtures."""
import gzip
import importlib.util
import json
import os
from pathlib import Path
import shutil
import subprocess
import tempfile
import unittest

ROOT = Path(__file__).resolve().parents[1]
SPEC = importlib.util.spec_from_file_location("indexed_plugin_cache", ROOT / "scripts/indexed_plugin_cache.py")
HELPER = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(HELPER)


@unittest.skipUnless(shutil.which("bgzip") and shutil.which("tabix"), "htslib tools required")
class IndexedPluginCache(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory(
            prefix="priva_indexed_plugin_", dir=os.environ.get("PRIVA_TEST_TMPDIR", os.environ.get("TMPDIR")))
        self.addCleanup(self.temporary.cleanup)
        self.folder = Path(self.temporary.name)

    def bgzf(self, content, name="table.tsv.gz"):
        path = self.folder / name
        with path.open("wb") as output:
            subprocess.run(["bgzip", "-c"], input=content, stdout=output, check=True)
        return path

    def index(self, path, *options):
        subprocess.run(["tabix", *options, str(path)], check=True)

    def test_mavedb_stale_zero_based_header_skip_and_boundary_queries_preserved(self):
        content = (b"unprefixed skipped header\n!publisher comment\n"
                   b"chr1\t0\t1\tA\tC\nchr1\t1\t2\tG\tT\nchr1\t2\t4\tAA\tC\n")
        path = self.bgzf(content)
        self.index(path, "-0", "-s", "1", "-b", "2", "-e", "3", "-c", "!", "-S", "1")
        index = Path(str(path) + ".tbi")
        settings = HELPER.index_settings(index)
        original = HELPER.stat_identity(path), path.read_bytes()
        queries = {pos: subprocess.check_output(["tabix", str(path), f"chr1:{pos}-{pos}"])
                   for pos in (1, 2, 3, 4, 5)}
        self.assertEqual(queries[1], b"chr1\t0\t1\tA\tC\n")
        self.assertEqual(queries[2], b"chr1\t1\t2\tG\tT\n")
        os.utime(index, (1, 1))
        result = HELPER.ensure_index(input_file=path, kind="mavedb")
        self.assertTrue(result["repaired"])
        self.assertEqual(HELPER.index_settings(index), settings)
        self.assertEqual(settings, (0x10000, 1, 2, 3, ord("!"), 1))
        self.assertEqual((HELPER.stat_identity(path), path.read_bytes()), original)
        for pos, expected in queries.items():
            self.assertEqual(subprocess.check_output(["tabix", str(path), f"chr1:{pos}-{pos}"]), expected)

    def test_missing_mavedb_index_uses_publisher_one_based_columns(self):
        path = self.bgzf(b"#header\nchr1\t1\t1\tA\tC\nchr1\t2\t2\tG\tT\n")
        self.assertTrue(HELPER.ensure_index(input_file=path, kind="mavedb")["repaired"])
        self.assertEqual(HELPER.index_settings(Path(str(path) + ".tbi")), (0, 1, 2, 3, 35, 0))
        self.assertEqual(subprocess.check_output(["tabix", str(path), "chr1:1-1"]), b"chr1\t1\t1\tA\tC\n")
        self.assertEqual(subprocess.check_output(["tabix", str(path), "chr1:2-2"]), b"chr1\t2\t2\tG\tT\n")

    def test_missing_index_does_not_guess_zero_based_custom_mavedb_coordinates(self):
        path = self.bgzf(b"chr1\t0\t1\tA\tC\n")
        with self.assertRaisesRegex(ValueError, "coordinates"):
            HELPER.ensure_index(input_file=path, kind="mavedb")
        self.assertFalse(Path(str(path) + ".tbi").exists())

    def test_current_index_and_table_are_untouched(self):
        path = self.bgzf(b"#chrom\tpos\nchr1\t1\nchr1\t2\n")
        HELPER.ensure_index(input_file=path, kind="primateai")
        index = Path(str(path) + ".tbi")
        before = HELPER.stat_identity(path), HELPER.stat_identity(index), index.read_bytes()
        self.assertFalse(HELPER.ensure_index(input_file=path, kind="primateai")["repaired"])
        self.assertEqual((HELPER.stat_identity(path), HELPER.stat_identity(index), index.read_bytes()), before)

    def test_corrupt_index_repaired_without_rewriting_supplied_table(self):
        path = self.bgzf(b"#chrom\tpos\nchr1\t1\nchr1\t2\n")
        Path(str(path) + ".tbi").write_bytes(b"corrupt")
        before = HELPER.stat_identity(path), path.read_bytes()
        self.assertTrue(HELPER.ensure_index(input_file=path, kind="splicevault")["repaired"])
        self.assertEqual((HELPER.stat_identity(path), path.read_bytes()), before)

    def test_gzip_empty_and_malformed_selected_inputs_fail(self):
        for content in (gzip.compress(b"chr1\t1\n"), b"", b"invalid"):
            path = self.folder / "invalid.gz"
            path.write_bytes(content)
            with self.assertRaises(ValueError):
                HELPER.ensure_index(input_file=path, kind="primateai")

    def raw_primate(self):
        content = b"".join(f"#publisher metadata {i}\n".encode() for i in range(11))
        content += (b"chr\tpos\tref\talt\tscore\nchr10\t2\tA\tT\t0.1\n\n"
                    b"chr2\t10\tG\tT\t0.5\nchr2\t2\tC\tA\t0.3\n"
                    b"chr1\t1\tA\tG\t0.7\nchr2\t2\tC\tG\t0.2\n")
        path = self.folder / "manual publisher.tsv.gz"
        with gzip.open(path, "wb") as stream:
            stream.write(content)
        return path

    def test_primate_matches_original_pipeline_and_binds_source_assembly(self):
        source = self.raw_primate()
        before = HELPER.stat_identity(source), source.read_bytes()
        reference = subprocess.check_output([
            "bash", "-o", "pipefail", "-c",
            '''gunzip -cf "$1" | sed '12s/.*/#&/' | sed '/^$/d' | awk 'NR<12{print $0;next}{print $0 | "sort -k1,1 -k 2,2n -V"}' ''',
            "test", str(source)])
        result = HELPER.format_primate(input_file=source, cache_dir=self.folder / "cache", assembly="GRCh38")
        output = Path(result["path"])
        self.assertEqual(gzip.decompress(output.read_bytes()), reference)
        receipt = json.loads((output.parent / "source.ready.json").read_text())
        self.assertEqual(receipt["source"]["assembly"], "GRCh38")
        self.assertEqual(receipt["source"]["sha256"], HELPER.sha256(source))
        self.assertEqual(receipt["source"]["path"], str(source))
        self.assertEqual((HELPER.stat_identity(source), source.read_bytes()), before)
        self.assertFalse(list(self.folder.glob("*.sorted*")))
        self.assertEqual(subprocess.check_output(["tabix", str(output), "chr2:2-2"]),
                         b"chr2\t2\tC\tA\t0.3\nchr2\t2\tC\tG\t0.2\n")

    def test_primate_same_source_reuses_and_changed_declared_assembly_is_distinct(self):
        source = self.raw_primate()
        args = dict(input_file=source, cache_dir=self.folder / "cache", assembly="hg38")
        first = Path(HELPER.format_primate(**args)["path"])
        before = HELPER.stat_identity(first)
        self.assertEqual(Path(HELPER.format_primate(**args)["path"]), first)
        self.assertEqual(HELPER.stat_identity(first), before)
        second = Path(HELPER.format_primate(**dict(args, assembly="hg19"))["path"])
        self.assertNotEqual(first, second)
        self.assertEqual(gzip.decompress(first.read_bytes()), gzip.decompress(second.read_bytes()))

    def test_primate_occupied_or_changed_derived_output_fails(self):
        source = self.raw_primate()
        args = dict(input_file=source, cache_dir=self.folder / "cache", assembly="hg38")
        output = Path(HELPER.format_primate(**args)["path"])
        output.write_bytes(b"changed derived output")
        with self.assertRaisesRegex(ValueError, "Conflicting"):
            HELPER.format_primate(**args)
        self.assertEqual(output.read_bytes(), b"changed derived output")


if __name__ == "__main__":
    unittest.main()
