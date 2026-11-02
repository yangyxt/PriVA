"""Exercise the actual estimator command with literal filesystem paths.

Only the bcftools boundary is a fixture. The source's shell, tr, wc and integer
conversion execute unchanged. Scientific VCF/table comparisons belong to the
real dependency test environment; these tests require only Python and POSIX tools.
"""
from __future__ import annotations

import ast
import json
import os
from pathlib import Path
import shlex
import subprocess
import sys
import tempfile
import unittest
from unittest import mock


SOURCE = Path(__file__).resolve().parents[1]/"scripts/combine_annotations.py"


class EstimatorPaths(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory()
        self.addCleanup(self.temporary.cleanup)
        self.root = Path(self.temporary.name)
        self.bin = self.root/"bin"
        self.bin.mkdir()
        self.record = self.root/"argv.json"
        executable = self.bin/"bcftools"
        executable.write_text("#!"+sys.executable+"\n"+'''import json, os, pathlib, sys
args = sys.argv[1:]
pathlib.Path(os.environ["ESTIMATOR_TEST_ARGV"]).write_text(json.dumps(args))
if len(args) != 4 or args[:3] != ["query", "-f", "%INFO/CSQ\\\\n"]:
    sys.exit(2)
sys.stdout.write(pathlib.Path(args[3]).read_text())
''')
        executable.chmod(0o755)
        self.environment = mock.patch.dict(os.environ, {
            "PATH": str(self.bin)+os.pathsep+os.environ.get("PATH", ""),
            "ESTIMATOR_TEST_ARGV": str(self.record),
        })
        self.environment.start()
        self.addCleanup(self.environment.stop)
        tree = ast.parse(SOURCE.read_text(), filename=str(SOURCE))
        function = next(node for node in tree.body
                        if isinstance(node, ast.FunctionDef) and node.name == "convert_vcf_to_tab")
        # Run the real allocation-estimate statements, stopping before the
        # unrelated dataframe/multiprocessing work and optional dependencies.
        statements = function.body[:2]
        self.assertEqual([node.targets[0].id for node in statements],
                         ["est_anno_count", "est_anno_count"])
        self.command = compile(ast.Module(body=statements, type_ignores=[]), str(SOURCE), "exec")

    def estimate(self, path):
        namespace = {"input_vcf": path, "subprocess": subprocess, "shlex": shlex}
        exec(self.command, namespace)
        return namespace["est_anno_count"]

    def test_literal_paths_preserve_estimated_count_and_the_complete_argument(self):
        # There is no invocation of the text in these filenames. In particular
        # command substitution must remain literal data passed to bcftools.
        names = ["input.vcf", "input space.vcf", "input\tname.vcf", "input\nname.vcf",
                 "input'quote.vcf", 'input"quote.vcf', "input;false.vcf",
                 "input$(false).vcf", "input`false`.vcf", "input*[x]?.vcf",
                 "input\\name.vcf", "inputé.vcf", "-input.vcf"]
        for name in names:
            with self.subTest(name=name):
                path = self.root/name
                path.write_text("a|b,c|d\nx|y\n")
                self.assertEqual(self.estimate(path), 3)
                self.assertEqual(json.loads(self.record.read_text()),
                                 ["query", "-f", "%INFO/CSQ\\n", str(path)])

    def test_existing_empty_and_unterminated_line_count_semantics(self):
        for text, expected in (("", 0), ("one", 0), ("one,two", 1),
                               (".\n", 1), ("a,b\n\n", 3)):
            with self.subTest(text=text):
                path = self.root/"input.vcf"
                path.write_text(text)
                self.assertEqual(self.estimate(path), expected)

    def test_missing_pipeline_tools_retain_the_integer_conversion_failure(self):
        with mock.patch.dict(os.environ, {"PATH": ""}):
            with self.assertRaisesRegex(ValueError, "invalid literal for int"):
                self.estimate(self.root/"input.vcf")


if __name__ == "__main__":
    unittest.main()
