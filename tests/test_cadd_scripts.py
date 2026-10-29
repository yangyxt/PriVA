"""PriVA's CADD prescore lookup and MMSplice copies: installation and equivalence.

The lookup fixture builds a small public-format prescored table; the reference
below repeats the upstream per-variant loop of CADD's extract_scored.py.
"""
import importlib.util
import os
from pathlib import Path
import random
import shutil
import subprocess
import sys
import tempfile
import unittest
from unittest import mock

try:
    import pysam
except ImportError:
    pysam = None


ROOT = Path(__file__).resolve().parents[1]
COPIES = ROOT / "scripts/cadd_scripts"
LOOKUP = COPIES / "extract_scored.py"


def load(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


installer = load("test_cadd_scripts_install_module", ROOT / "scripts/cadd_scripts_install.py")
index_builder = load("test_cadd_prescore_index_module", ROOT / "scripts/cadd_prescore_index.py")


class InstallCopies(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory(prefix="cadd_scripts_", dir=os.environ.get("PRIVA_TEST_TMPDIR"))
        self.addCleanup(self.temp.cleanup)
        self.source = Path(self.temp.name) / "CADD-scripts"
        (self.source / "src/scripts/lib/tools").mkdir(parents=True)
        (self.source / "CADD.sh").write_text("#!/bin/sh\n")
        # Fixture upstream files stand in for the reviewed CADD bytes.
        table = {}
        for relative, (copy, _) in installer.FILES.items():
            payload = ("# fixture upstream " + relative + "\n").encode()
            path = self.source / relative
            path.write_bytes(payload)
            path.chmod(0o755 if relative.endswith("extract_scored.py") else 0o644)
            table[relative] = (copy, installer.blob_id(payload))
        patcher = mock.patch.dict(installer.FILES, table)
        patcher.start()
        self.addCleanup(patcher.stop)

    def states(self):
        return {relative: installer.state(self.source, relative) for relative in installer.FILES}

    def test_install_replaces_upstream_bytes_and_keeps_mode(self):
        installer.install(self.source)
        self.assertEqual(set(self.states().values()), {"priva"})
        for relative, (copy, _) in installer.FILES.items():
            path = self.source / relative
            self.assertEqual(path.read_bytes(), (COPIES / copy).read_bytes())
            self.assertEqual(path.stat().st_mode & 0o777, 0o755 if relative.endswith("extract_scored.py") else 0o644)
        installer.install(self.source)
        self.assertEqual(set(self.states().values()), {"priva"})

    def test_earlier_priva_copy_is_updated(self):
        relative = "src/scripts/extract_scored.py"
        earlier = b"# an earlier PriVA copy\n"
        (self.source / relative).write_bytes(earlier)
        self.assertEqual(self.states()[relative], "other")
        with mock.patch.object(installer, "PREVIOUS", {installer.blob_id(earlier)}):
            self.assertEqual(self.states()[relative], "previous priva")
            installer.install(self.source)
        self.assertEqual(self.states()[relative], "priva")

    def test_other_content_is_left_unchanged(self):
        relative = "src/scripts/extract_scored.py"
        (self.source / relative).write_bytes(b"# a local edit\n")
        installer.install(self.source)
        self.assertEqual((self.source / relative).read_bytes(), b"# a local edit\n")
        self.assertEqual(self.states()[relative], "other")
        self.assertEqual(self.states()["src/scripts/lib/tools/MMSplice.py"], "priva")

    def test_command_reports_state_and_refuses_a_folder_without_cadd_sh(self):
        command = [sys.executable, str(ROOT / "scripts/cadd_scripts_install.py"), "--source"]
        result = subprocess.run(command + [str(self.source), "--state"], capture_output=True, text=True)
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertEqual(len(result.stdout.splitlines()), 2)
        result = subprocess.run(command + [str(self.source / "src")], capture_output=True, text=True)
        self.assertNotEqual(result.returncode, 0)
        self.assertIn("CADD.sh was not found", result.stderr)


class PrivaCopies(unittest.TestCase):
    def test_copies_compile_and_carry_the_cadd_notice(self):
        for name in ("extract_scored.py", "MMSplice.py"):
            path = COPIES / name
            compile(path.read_text(), str(path), "exec")
            self.assertIn("Hudson-Alpha Institute for", path.read_text())
        self.assertTrue((COPIES / "CADD-LICENSE").is_file())

    def test_lookup_is_ascii_for_upstream_python2(self):
        LOOKUP.read_bytes().decode("ascii")


def reference_lookup(table, lines):
    """Upstream extract_scored.py's loop: one tabix query per variant."""
    tabix = pysam.TabixFile(str(table))
    novel, found, messages = [], [], 0
    fpos, fref, falt = 1, 2, 3
    for header in tabix.header:
        found.append(header + "\n")
        fields = header.split("\t")
        if "Ref" in fields:
            fref = fields.index("Ref")
            if "Alt" in fields:
                falt = fields.index("Alt")
    for line in lines:
        line = line.rstrip("\n\r")
        if line.startswith("#"):
            novel.append(line + "\n")
            continue
        try:
            fields = line.split("\t")
            hit = False
            pos = int(fields[1])
            for row in tabix.fetch(fields[0], pos - 1, pos):
                values = row.rstrip().split("\t")
                if values[fref] == fields[-2] and values[falt] == fields[-1] and values[fpos] == fields[1]:
                    found.append(row + "\n")
                    hit = True
            if not hit:
                novel.append(line + "\n")
        except ValueError:
            messages += 1
            novel.append(line + "\n")
    return "".join(novel), "".join(found), messages


@unittest.skipUnless(pysam, "pysam is required for the prescored table fixture")
class LookupEquivalence(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.temp = tempfile.TemporaryDirectory(prefix="cadd_lookup_", dir=os.environ.get("PRIVA_TEST_TMPDIR"))
        root = Path(cls.temp.name)
        generator = random.Random(20261028)
        rows = []
        # Dense stretches, isolated positions further apart than SPAN_GAP and
        # repeated rows, on two chromosomes in tabix order.
        for chrom, starts in (("1", (1000, 5000, 400000, 2000000)), ("2", (50, 900000))):
            for start in starts:
                for pos in range(start, start + 1200, generator.choice((1, 2, 7))):
                    for alt in "ACGT"[: generator.randint(1, 3)]:
                        rows.append((chrom, pos, "N", alt, generator.random()))
            rows.append((chrom, starts[-1] + 10, "N", "A", 0.5))
        rows.sort(key=lambda row: (row[0], row[1]))
        plain = root / "prescored.tsv"
        with plain.open("w") as stream:
            stream.write("## CADD fixture\n#Chrom\tPos\tRef\tAlt\tRawScore\n")
            for chrom, pos, ref, alt, score in rows:
                stream.write("%s\t%d\t%s\t%s\t%.4f\n" % (chrom, pos, ref, alt, score))
        cls.table = Path(pysam.tabix_index(str(plain), seq_col=0, start_col=1, end_col=1,
                                           meta_char="#", force=True))
        lines = ["##fileformat=VCFv4.2", "#CHROM\tPOS\tID\tREF\tALT"]
        for chrom, pos, _, alt, _ in generator.sample(rows, min(len(rows), 2500)):
            lines.append("%s\t%d\t.\tN\t%s" % (chrom, pos, alt))
            lines.append("%s\t%d\t.\tN\tTT" % (chrom, pos))
        lines += ["1\t1001\t.\tN\tA", "1\t1001\t.\tN\tA", "1\t01001\t.\tN\tA", "1\t0\t.\tN\tA",
                  "1\t-5\t.\tN\tA", "1\tseven\t.\tN\tA", "Z\t1001\t.\tN\tA", "2\t999999999\t.\tN\tA",
                  "#a header in the middle", "2\t60\t.\tN\tC"]
        generator.shuffle(lines)
        cls.lines = [line + "\n" for line in lines]

    @classmethod
    def tearDownClass(cls):
        cls.temp.cleanup()

    def run_lookup(self, found_file=True, interpreter=sys.executable, table=None):
        root = Path(self.temp.name)
        found = root / "found.tsv"
        found.unlink(missing_ok=True)
        command = [interpreter, str(LOOKUP), "--header", "-p", str(table or self.table)]
        if found_file:
            command.append("--found_out=" + str(found))
        result = subprocess.run(command, input="".join(self.lines), capture_output=True, text=True, timeout=120)
        self.assertEqual(result.returncode, 0, result.stderr)
        return result.stdout, found.read_text() if found_file else None, result.stderr.count("Encountered uncovered chromosome")

    def test_rows_and_order_match_upstream(self):
        # CADD runs the lookup under its own Python 2.7 environment when given.
        interpreters = [sys.executable] + [os.environ[key] for key in ("PRIVA_TEST_CADD_PYTHON2",) if os.environ.get(key)]
        expected = reference_lookup(self.table, self.lines)
        for interpreter in interpreters:
            with self.subTest(interpreter=interpreter):
                novel, found, messages = self.run_lookup(interpreter=interpreter)
                self.assertEqual((novel, found, messages), expected)
                self.assertGreater(found.count("\n"), 1000)
                self.assertGreater(messages, 0)

    def test_combined_output_keeps_upstream_interleaving(self):
        # Without --found_out upstream writes found rows and novel lines to stdout in input order.
        stdout, _, _ = self.run_lookup(found_file=False)
        tabix = pysam.TabixFile(str(self.table))
        expected = "".join(header + "\n" for header in tabix.header)
        for line in self.lines:
            single_novel, single_found, _ = reference_lookup(self.table, [line])
            expected += single_found[len("".join(header + "\n" for header in tabix.header)):] + single_novel
        self.assertEqual(stdout, expected)

    def test_batches_split_a_stretch_without_changing_rows(self):
        root = Path(self.temp.name)
        small = root / "extract_scored_small_batch.py"
        small.write_text(LOOKUP.read_text().replace("BATCH_SIZE = 20000", "BATCH_SIZE = 7")
                         .replace("SPAN_GAP = 256 if position_index else 4096", "SPAN_GAP = 100000"))
        found = root / "found_small.tsv"
        result = subprocess.run([sys.executable, str(small), "--header", "-p", str(self.table),
                                 "--found_out=" + str(found)], input="".join(self.lines),
                                capture_output=True, text=True, timeout=120)
        self.assertEqual(result.returncode, 0, result.stderr)
        novel, rows, messages = reference_lookup(self.table, self.lines)
        self.assertEqual((result.stdout, found.read_text()), (novel, rows))

    def test_malformed_line_stops_after_writing_earlier_results(self):
        lines = ["1\t1001\t.\tN\tA\n", "only-one-field\n", "1\t1003\t.\tN\tA\n"]
        root = Path(self.temp.name)
        found = root / "found_stop.tsv"
        result = subprocess.run([sys.executable, str(LOOKUP), "-p", str(self.table), "--found_out=" + str(found)],
                                input="".join(lines), capture_output=True, text=True, timeout=60)
        self.assertNotEqual(result.returncode, 0)
        self.assertIn("IndexError", result.stderr)
        self.assertEqual(found.read_text(), reference_lookup(self.table, lines[:1])[1].split("#Chrom\tPos\tRef\tAlt\tRawScore\n")[-1])


@unittest.skipUnless(pysam, "pysam is required for the prescored table fixture")
class PositionIndexLookup(unittest.TestCase):
    """A multi-block table indexed by scripts/cadd_prescore_index.py."""

    @classmethod
    def setUpClass(cls):
        cls.temp = tempfile.TemporaryDirectory(prefix="cadd_posidx_", dir=os.environ.get("PRIVA_TEST_TMPDIR"))
        root = Path(cls.temp.name)
        generator = random.Random(20261029)
        rows = []
        # About 250 BGZF blocks: dense runs, sparse stretches and repeated rows.
        for chrom in ("1", "2", "X"):
            position = 10
            while position < 3000000:
                for alt in "ACG"[: generator.randint(1, 3)]:
                    rows.append((chrom, position, alt))
                position += generator.choice((1, 1, 1, 2, 3, 5, 40, 900))
            rows.append((chrom, position - 1, "T"))
        rows.sort(key=lambda row: (row[0], row[1]))
        plain = root / "prescored.tsv"
        padding = "x" * 120
        with plain.open("w") as stream:
            stream.write("## CADD fixture\n#Chrom\tPos\tRef\tAlt\tAnno\n")
            for chrom, position, alt in rows:
                stream.write("%s\t%d\tN\t%s\t%s%d\n" % (chrom, position, alt, padding, generator.randint(0, 9)))
        cls.table = Path(pysam.tabix_index(str(plain), seq_col=0, start_col=1, end_col=1,
                                           meta_char="#", force=True))
        built = subprocess.run([sys.executable, str(ROOT / "scripts/cadd_prescore_index.py"), "build",
                                str(cls.table), "--threads", "3", "--every", "1"],
                               capture_output=True, text=True, timeout=300, check=True)
        cls.index = Path(str(cls.table) + ".posidx")
        cls.count = int(built.stdout.split(" with ")[1].split()[0])
        lines = []
        for chrom, position, alt in generator.sample(rows, min(len(rows), 4000)):
            lines.append("%s\t%d\t.\tN\t%s" % (chrom, position, alt))
            lines.append("%s\t%d\t.\tN\tTT" % (chrom, position + 1))
        for chrom in ("1", "2", "X"):
            # Chromosome starts precede the first indexed row and use tabix.
            lines += ["%s\t10\t.\tN\tA" % chrom, "%s\t11\t.\tN\tA" % chrom, "%s\t2999999\t.\tN\tT" % chrom]
        lines.sort(key=lambda line: (line.split("\t")[0], int(line.split("\t")[1])))
        cls.lines = [line + "\n" for line in lines]

    @classmethod
    def tearDownClass(cls):
        cls.temp.cleanup()

    def lookup(self, interpreter=sys.executable, table=None):
        table = table or self.table
        found = Path(self.temp.name) / "found.tsv"
        found.unlink(missing_ok=True)
        result = subprocess.run([interpreter, str(LOOKUP), "--header", "-p", str(table), "--found_out=" + str(found)],
                                input="".join(self.lines), capture_output=True, text=True, timeout=300)
        self.assertEqual(result.returncode, 0, result.stderr)
        return result.stdout, found.read_text(), result.stderr

    def test_index_covers_many_blocks_and_matches_its_table(self):
        self.assertGreater(self.count, 150)
        self.assertTrue(index_builder.check(self.table))

    def test_rows_and_order_match_upstream_with_the_index(self):
        expected = reference_lookup(self.table, self.lines)
        interpreters = [sys.executable] + [os.environ[key] for key in ("PRIVA_TEST_CADD_PYTHON2",) if os.environ.get(key)]
        for interpreter in interpreters:
            with self.subTest(interpreter=interpreter):
                novel, found, errors = self.lookup(interpreter)
                self.assertNotIn("Ignoring", errors)
                self.assertEqual((novel, found), expected[:2])
                self.assertGreater(found.count("\n"), 4000)

    def test_index_of_another_table_is_ignored(self):
        other = Path(self.temp.name) / "other" / self.table.name
        other.parent.mkdir(exist_ok=True)
        shutil.copyfile(self.table, other)
        shutil.copyfile(str(self.table) + ".tbi", str(other) + ".tbi")
        # The same index bytes, but recording a different table size.
        payload = bytearray(Path(str(self.index)).read_bytes())
        start = len(index_builder.MAGIC)
        payload[start:start + 8] = (other.stat().st_size + 1).to_bytes(8, "little")
        Path(str(other) + ".posidx").write_bytes(bytes(payload))
        self.assertFalse(index_builder.check(other))
        novel, found, errors = self.lookup(table=other)
        self.assertIn("Ignoring", errors)
        self.assertEqual((novel, found), reference_lookup(self.table, self.lines)[:2])


if __name__ == "__main__":
    unittest.main()
