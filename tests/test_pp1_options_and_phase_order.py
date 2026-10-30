from __future__ import annotations

import ast
from pathlib import Path
import sys
import unittest
from unittest import mock

import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
SCRIPTS = ROOT / "scripts"
sys.path.insert(0, str(SCRIPTS))

import determine_phase  # noqa: E402


def _main_block(tree: ast.Module) -> ast.If:
    for node in tree.body:
        if isinstance(node, ast.If) and "__main__" in ast.unparse(node.test):
            return node
    raise AssertionError("acmg_criteria_assign.py has no __main__ block")


class Pp1CommandLineTest(unittest.TestCase):
    def test_acmg_command_line_accepts_and_forwards_the_pp1_inputs(self) -> None:
        shell = (SCRIPTS / "prioritization_vcf_per_fam.sh").read_text(encoding="utf-8")
        self.assertIn('--pp1_vcf ${pp1_vcf} --pp1_ped ${pp1_ped}', shell)

        tree = ast.parse((SCRIPTS / "acmg_criteria_assign.py").read_text(encoding="utf-8"))
        main = _main_block(tree)
        options = {
            call.args[0].value
            for call in ast.walk(main)
            if isinstance(call, ast.Call)
            and ast.unparse(call.func) == "parser.add_argument"
            and call.args
            and isinstance(call.args[0], ast.Constant)
        }
        self.assertLessEqual({"--pp1_vcf", "--pp1_ped"}, options)

        (assign_call,) = [
            call
            for call in ast.walk(main)
            if isinstance(call, ast.Call) and ast.unparse(call.func) == "ACMG_criteria_assign"
        ]
        forwarded = {keyword.arg: ast.unparse(keyword.value) for keyword in assign_call.keywords}
        self.assertEqual(forwarded.get("pp1_vcf"), "args.pp1_vcf")
        self.assertEqual(forwarded.get("pp1_ped"), "args.pp1_ped")


class _PoolReturningGenesInReverse:
    """Stands in for multiprocessing.Pool: imap keeps the order, imap_unordered reverses it."""

    def __init__(self, processes: int) -> None:
        self.processes = processes

    def imap(self, function, iterable, chunksize=1):
        return [function(item) for item in iterable]

    def imap_unordered(self, function, iterable, chunksize=1):
        return list(reversed([function(item) for item in iterable]))

    def close(self) -> None:
        pass

    def join(self) -> None:
        pass


class PhaseMergeOrderTest(unittest.TestCase):
    def test_phase_of_a_variant_shared_by_two_genes_does_not_depend_on_process_count(self) -> None:
        variants = pd.DataFrame(
            {
                "variant_id": ["chr1:100:A-G", "chr1:100:A-G", "chr1:200:C-T"],
                "Gene": ["GENE_A", "GENE_B", "GENE_B"],
                "P1": ["0/1", "0/1", "0/1"],
            }
        )
        phase_by_gene = {"GENE_A": "cis", "GENE_B": "trans"}

        def phases(processes: int) -> list:
            with mock.patch.object(determine_phase, "build_patient_info", lambda pedigree, df: {"P1": {}}), \
                    mock.patch.object(determine_phase, "_can_resolve_phase", lambda df, info: True), \
                    mock.patch.object(determine_phase, "process_gene_wrapper",
                                      lambda item: {"P1": {"chr1:100:A-G": phase_by_gene[item[0]]}}), \
                    mock.patch.object(determine_phase.mp, "Pool", _PoolReturningGenesInReverse):
                _, _, table = determine_phase.determine_cis_trans_relationships(
                    variants.copy(), [False, False, True], pd.DataFrame(), num_processes=processes
                )
            return table["P1_phase_patho"].tolist()

        one_process = phases(1)
        self.assertEqual(one_process[:2], ["trans", "trans"])
        self.assertEqual(phases(4), one_process)


if __name__ == "__main__":
    unittest.main()
