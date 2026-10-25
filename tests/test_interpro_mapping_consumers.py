"""Tiny real producer/reader tests, with explicitly isolated routing boundaries.

The real InterPro XML/VCF producer writes the synthetic fixture. Its bytes and
an equivalent gzip stream reach the actual mapping-reader statements extracted
from DAS, PVS1 and splicing. AST extraction avoids running unrelated scientific
stages; it does not replace deserialization. CLI/PM1 forwarding uses labelled
capture boundaries, not a scientific DAS build or whole-installation acceptance.
"""
import argparse
import ast
import gzip
import importlib.util
import json
import logging
import os
from pathlib import Path
import pickle
import subprocess
import sys
import tempfile
import types
import unittest
from unittest import mock


ROOT = Path(__file__).resolve().parents[1]
SCRIPTS = ROOT / "scripts"
SPEC = importlib.util.spec_from_file_location("interpro_consumer_io", SCRIPTS / "interpro_mapping_io.py")
IO = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(IO)
CONSUMERS = {
    "am_pick_intolerant_domains.py": ("interpro_map_pickle", "interpro_map_dict"),
    "acmg_pvs1_null_variant.py": ("interpro_entry_map_dict",),
    "splicing_var_analysis.py": ("interpro_entry_map_dict",),
}
XML = b'''<interprodb><release version="synthetic"/>
<interpro id="IPR000001" type="Domain" short_name="First"><name>First name</name>
<class_list><classification class_type="GO" id="GO:0000001"><category>category</category><description>description</description></classification></class_list>
<abstract>First <p>abstract</p> text.</abstract><member_list>
<db_xref db="Pfam" dbkey="PF00001"/><db_xref db="Pfam" dbkey="PF00001"/>
<db_xref db="SMART" dbkey="SM00001"/>
</member_list></interpro>
<interpro id="IPR000002"><name>Second name</name><member_list><db_xref db="PFAM" dbkey="PF00001"/></member_list></interpro>
</interprodb>'''


def source_tree(name):
    path = SCRIPTS / name
    return ast.parse(path.read_text(), filename=str(path))


def execute_nodes(nodes, namespace, filename):
    """Execute exact source nodes only; callers explicitly label boundaries."""
    code = compile(ast.fix_missing_locations(ast.Module(body=nodes, type_ignores=[])), filename, "exec")
    exec(code, namespace)


class InterProMappingConsumers(unittest.TestCase):
    def setUp(self):
        retained = os.environ.get("PRIVA_INTERPRO_CONSUMER_FIXTURE_ROOT")
        if retained:
            self.folder = Path(retained) / self._testMethodName
            self.folder.mkdir(parents=True, exist_ok=False)
        else:
            temporary = tempfile.TemporaryDirectory(prefix="interpro_consumer_fixture_", dir=os.environ.get(
                "PRIVA_TEST_TMPDIR", os.environ.get("TMPDIR")))
            self.addCleanup(temporary.cleanup)
            self.folder = Path(temporary.name)
        # Deliberately include ordered keys and duplicate values.
        self.mapping = {("PFAM", "PF00001"): [("IPR1", "Domain"), ("IPR1", "Domain"),
                                              ("IPR2", "Repeat")],
                        ("SMART", "SM00001"): [("IPR3", "Domain")]}
        self.raw = pickle.dumps(self.mapping, protocol=pickle.HIGHEST_PROTOCOL)
        self.fake_scripts = self.folder / "standalone_source/scripts"
        self.fake_scripts.mkdir(parents=True)
        self.bundle = self.fake_scripts.parent / "data/InterPro/Interpro_entry_mapping.pkl.gz"
        self.bundle.parent.mkdir(parents=True)
        self.bundle.write_bytes(gzip.compress(self.raw, mtime=0))

    def snapshot(self, path):
        value = path.stat()
        return (path.read_bytes(), value.st_dev, value.st_ino, value.st_size,
                value.st_mtime_ns, value.st_ctime_ns, value.st_mode, value.st_nlink)

    def write(self, name, content):
        path = self.folder / name
        path.write_bytes(content)
        return path

    def consumer_read(self, name, path):
        """Real deserialization at an isolated actual source-reader boundary."""
        tree = source_tree(name)
        imports = [node for node in tree.body if isinstance(node, ast.ImportFrom)
                   and node.module == "interpro_mapping_io"]
        self.assertEqual(len(imports), 1)
        self.assertIn("load_interpro_mapping", [alias.name for alias in imports[0].names])
        nodes = []
        for target in CONSUMERS[name]:
            matches = [node for node in ast.walk(tree) if isinstance(node, ast.Assign)
                       and any(isinstance(item, ast.Name) and item.id == target for item in node.targets)]
            self.assertEqual(len(matches), 1, (name, target))
            nodes.append(matches[0])
        namespace = {"os": os, "gzip": gzip, "pickle": pickle,
                     "load_interpro_mapping": IO.load_interpro_mapping,
                     "interpro_mapping_pickle": None if path is None else str(path),
                     "interpro_entry_map_pkl": None if path is None else str(path),
                     "self_directory": str(self.fake_scripts)}
        execute_nodes(nodes, namespace, str(SCRIPTS / name))
        return namespace[CONSUMERS[name][-1]]

    def all_readers(self, path, expected):
        before = self.snapshot(path)
        self.assertEqual(IO.load_interpro_mapping(path), expected)
        for name in CONSUMERS:
            with self.subTest(consumer=name):
                value = self.consumer_read(name, path)
                self.assertEqual(value, expected)
                self.assertEqual(list(value), list(expected))
                for key in expected:
                    self.assertEqual(value[key], expected[key])
        self.assertEqual(self.snapshot(path), before)

    def producer(self):
        xml = self.write("interpro.xml.gz", gzip.compress(XML, mtime=0))
        vcf = self.folder / "synthetic.vcf"
        vcf.write_text('##fileformat=VCFv4.2\n##contig=<ID=chr1,length=100>\n'
                       '##INFO=<ID=CSQ,Number=.,Type=String,Description="Synthetic fixture. Format: Allele|DOMAINS">\n'
                       '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n'
                       'chr1\t1\t.\tA\tC\t.\t.\tCSQ=C|Pfam:PF00001\n')
        output = self.folder / "mapping.pkl.gz"
        command = [sys.executable, str(SCRIPTS / "protein_domain_mapping.py"),
                   "--xml", str(xml), "--vcf", str(vcf), "--mapping_output", str(output),
                   "--report", str(self.folder / "comparison.txt")]
        result = subprocess.run(command, text=True, capture_output=True, timeout=30)
        (self.folder / "producer.command.json").write_text(json.dumps(command, indent=2) + "\n")
        (self.folder / "producer.stdout.log").write_text(result.stdout)
        (self.folder / "producer.stderr.log").write_text(result.stderr)
        self.assertEqual(result.returncode, 0, result.stderr)
        return output

    def assert_corrupt_unchanged(self, content, suffix=".pkl.gz"):
        path = self.write("selected" + suffix, content)
        before = self.snapshot(path)
        with self.assertRaises((OSError, EOFError, ValueError, pickle.PickleError)):
            IO.load_interpro_mapping(path)
        self.assertEqual(self.snapshot(path), before)

    def assert_ordinary_gzip_read_parity(self, content):
        """Do not add trailer/content validation after the ordinary pickle read."""
        path = self.write("selected.pkl.gz", content)
        before = self.snapshot(path)
        try:
            with gzip.open(path, "rb") as handle:
                expected = pickle.load(handle)
        except (OSError, EOFError, ValueError, pickle.PickleError) as error:
            with self.assertRaises(type(error)):
                IO.load_interpro_mapping(path)
        else:
            self.assertEqual(IO.load_interpro_mapping(path), expected)
        self.assertEqual(self.snapshot(path), before)

    def test_real_producer_raw_gz_suffix_reaches_all_actual_readers(self):
        path = self.producer()
        self.assertEqual(path.read_bytes()[:2], bytes((0x80, pickle.HIGHEST_PROTOCOL)))
        with path.open("rb") as handle:
            expected = pickle.load(handle)
        self.assertEqual(list(expected), [("PFAM", "PF00001"), ("SMART", "SM00001")])
        self.assertEqual([item[0] for item in expected[("PFAM", "PF00001")]],
                         ["IPR000001", "IPR000001", "IPR000002"])
        self.all_readers(path, expected)

    def test_real_producer_gzip_equivalent_reaches_all_actual_readers(self):
        path = self.producer()
        raw = path.read_bytes()
        equivalent = self.write("gzip-equivalent.without-gz-suffix", gzip.compress(raw, mtime=0))
        self.all_readers(equivalent, pickle.loads(raw))

    def test_raw_bytes_ignore_filename_suffix(self):
        for name in ("plain.pkl", "plain.pkl.gz", "plain.data"):
            with self.subTest(name=name):
                self.all_readers(self.write(name, self.raw), self.mapping)

    def test_gzip_bytes_ignore_filename_suffix(self):
        for name in ("compressed.pkl", "compressed.pkl.gz", "compressed.data"):
            with self.subTest(name=name):
                self.all_readers(self.write(name, gzip.compress(self.raw, mtime=0)), self.mapping)

    def test_empty_selected_mapping_raises_unchanged(self):
        self.assert_corrupt_unchanged(b"")

    def test_invalid_plain_pickle_raises_unchanged(self):
        self.assert_corrupt_unchanged(b"not a pickle")

    def test_invalid_gzip_magic_stream_raises_without_raw_fallback(self):
        self.assert_corrupt_unchanged(b"\x1f\x8bnot a gzip stream")

    def test_gzip_containing_invalid_pickle_raises_unchanged(self):
        self.assert_corrupt_unchanged(gzip.compress(b"not a pickle", mtime=0))

    def test_truncated_gzip_footer_matches_ordinary_pickle_read(self):
        self.assert_ordinary_gzip_read_parity(gzip.compress(self.raw, mtime=0)[:-4])

    def test_invalid_gzip_checksum_matches_ordinary_pickle_read(self):
        data = bytearray(gzip.compress(self.raw, mtime=0))
        data[-8] ^= 1
        self.assert_ordinary_gzip_read_parity(bytes(data))

    def test_truncated_plain_pickle_raises_unchanged(self):
        self.assert_corrupt_unchanged(self.raw[:-1])

    def test_gzip_containing_truncated_pickle_raises_unchanged(self):
        self.assert_corrupt_unchanged(gzip.compress(self.raw[:-1], mtime=0))

    def test_missing_selected_mapping_does_not_guess_another_path(self):
        selected = self.folder / "missing.pkl.gz"
        with self.assertRaises(FileNotFoundError):
            IO.load_interpro_mapping(selected)
        self.assertFalse(selected.exists())

    def test_directory_selected_mapping_raises(self):
        with self.assertRaises(IsADirectoryError):
            IO.load_interpro_mapping(self.folder)

    def test_all_actual_reader_boundaries_propagate_corruption(self):
        selected = self.write("selected.pkl.gz", b"not a pickle")
        before = self.snapshot(selected)
        for name in CONSUMERS:
            with self.subTest(consumer=name):
                with self.assertRaises(pickle.UnpicklingError):
                    self.consumer_read(name, selected)
        self.assertEqual(self.snapshot(selected), before)

    def test_das_omitted_path_retains_legacy_bundled_default(self):
        before = self.snapshot(self.bundle)
        self.assertEqual(self.consumer_read("am_pick_intolerant_domains.py", None), self.mapping)
        self.assertEqual(self.snapshot(self.bundle), before)

    def test_das_explicit_invalid_path_never_falls_back_to_existing_bundle(self):
        before = self.snapshot(self.bundle)
        for selected in (self.folder / "missing.pkl.gz", ""):
            with self.subTest(selected=str(selected)):
                with self.assertRaises((FileNotFoundError, IsADirectoryError)):
                    self.consumer_read("am_pick_intolerant_domains.py", selected)
        self.assertEqual(self.snapshot(self.bundle), before)

    def cli_forward(self, arguments):
        # Labelled CLI capture boundary: real argparse/forwarding, no DAS build.
        tree = source_tree("am_pick_intolerant_domains.py")
        mains = [node for node in tree.body if isinstance(node, ast.If)
                 and ast.dump(node.test) == ast.dump(ast.parse("__name__ == '__main__'", mode="eval").body)]
        self.assertEqual(len(mains), 1)
        calls = []
        def completed_boundary(*args, **kwargs):
            calls.append((args, kwargs))
            return True

        namespace = {"ap": argparse, "analyze_domain_data": completed_boundary}
        argv = ["am_pick_intolerant_domains.py", "--pickle_file", "synthetic-scores.pkl",
                "--output_dir", str(self.folder)] + arguments
        with mock.patch.object(sys, "argv", argv):
            execute_nodes(mains[0].body, namespace, str(SCRIPTS / "am_pick_intolerant_domains.py"))
        self.assertEqual(len(calls), 1)
        return calls[0]

    def test_das_cli_named_option_forwards_selected_path(self):
        selected = str(self.folder / "selected mapping.pkl.gz")
        args, kwargs = self.cli_forward(["--interpro_mapping_pickle", selected])
        self.assertEqual(kwargs["interpro_mapping_pickle"], selected)
        self.assertEqual(args, ("synthetic-scores.pkl", str(self.folder), 62))

    def test_das_cli_hyphen_alias_forwards_selected_path(self):
        selected = str(self.folder / "selected.pkl.gz")
        self.assertEqual(self.cli_forward(["--interpro-mapping-pickle", selected])[1]["interpro_mapping_pickle"], selected)

    def test_das_cli_omitted_option_forwards_none(self):
        self.assertIsNone(self.cli_forward([])[1]["interpro_mapping_pickle"])

    def pm1_forward(self, config):
        # Labelled PM1 producer boundary; only path propagation/reload is tested.
        tree = source_tree("prepare_pm1_regions.py")
        functions = [node for node in tree.body if isinstance(node, ast.FunctionDef) and node.name == "run_das_pipeline"]
        self.assertEqual(len(functions), 1)
        calls = []

        def captured_analyze(*args, **kwargs):
            calls.append((args, kwargs))
            with (self.folder / "all_intolerant_domains.hg38.pkl").open("wb") as handle:
                pickle.dump({"synthetic-domain"}, handle)
            return True

        boundary = types.ModuleType("am_pick_intolerant_domains")
        boundary.analyze_domain_data = captured_analyze
        namespace = {"os": os, "pickle": pickle, "sys": sys,
                     "logger": logging.getLogger("interpro_consumer_fixture"),
                     "PRIVA_DIR": str(self.folder / "standalone_source"), "SCRIPT_DIR": str(SCRIPTS)}
        execute_nodes(functions, namespace, str(SCRIPTS / "prepare_pm1_regions.py"))
        with mock.patch.dict(sys.modules, {"am_pick_intolerant_domains": boundary}):
            result = namespace["run_das_pipeline"](config, str(self.folder), "hg38", 1)
        self.assertEqual(result, {"synthetic-domain"})
        self.assertEqual(len(calls), 1)
        return calls[0]

    def test_pm1_direct_das_caller_forwards_configured_mapping(self):
        selected = str(self.folder / "selected.pkl.gz")
        self.assertEqual(self.pm1_forward({"interpro_mapping_pickle": selected})[1]["interpro_mapping_pickle"], selected)

    def test_pm1_direct_das_caller_retains_omitted_default(self):
        self.assertIsNone(self.pm1_forward({})[1]["interpro_mapping_pickle"])

    def test_helper_cli_prints_only_metadata_without_rewriting(self):
        selected = self.write("selected.pkl.gz", self.raw)
        before = self.snapshot(selected)
        result = subprocess.run([sys.executable, str(SCRIPTS / "interpro_mapping_io.py"),
                                 "inspect", "--input", str(selected)], capture_output=True, text=True, timeout=10)
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertEqual(json.loads(result.stdout), {"input": str(selected), "value_type": "dict", "entries": 2})
        self.assertNotIn("PF00001", result.stdout)
        self.assertEqual(self.snapshot(selected), before)

    def test_helper_cli_corrupt_selected_input_fails_nonzero(self):
        selected = self.write("selected.pkl.gz", b"not a pickle")
        before = self.snapshot(selected)
        result = subprocess.run([sys.executable, str(SCRIPTS / "interpro_mapping_io.py"),
                                 "inspect", "--input", str(selected)], capture_output=True, text=True, timeout=10)
        self.assertNotEqual(result.returncode, 0)
        self.assertEqual(result.stdout, "")
        self.assertIn("ERROR:", result.stderr)
        self.assertEqual(self.snapshot(selected), before)


if __name__ == "__main__":
    unittest.main()
