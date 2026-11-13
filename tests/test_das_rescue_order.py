"""Protect domain-to-seed assignment before either DAS rescue dispatches work.

The actual candidate-selection functions are loaded without importing the large
scientific environment. The pool stops at task submission; these checks do not
stand in for numerical dip-test or complete-cache comparisons.
"""
import ast
import json
import os
from pathlib import Path
import statistics
import subprocess
import sys
from types import SimpleNamespace
import unittest


SOURCE = Path(__file__).resolve().parents[1] / "scripts/am_pick_intolerant_domains.py"
sys.path.insert(0, str(SOURCE.parents[1] / "rust/tests"))
from das_sorted_reference import render_source
FUNCTIONS = ("rescue_bimodal_go_mf_domains", "rescue_non_gomf_domains")


class Submitted(Exception):
    def __init__(self, tasks):
        self.tasks = tasks


class CapturePool:
    def __enter__(self):
        return self

    def __exit__(self, *_):
        return False

    def imap_unordered(self, _worker, tasks):
        raise Submitted(tasks)


def procedures():
    tree = ast.parse(render_source(SOURCE.read_bytes()))
    selected = [node for node in tree.body
                if isinstance(node, ast.FunctionDef) and node.name in FUNCTIONS]
    assert len(selected) == 2
    namespace = dict(np=SimpleNamespace(median=statistics.median),
                     logger=SimpleNamespace(info=lambda *_: None),
                     mp=SimpleNamespace(get_context=lambda method:
                         SimpleNamespace(Pool=lambda workers: CapturePool())),
                     _dip_rescue_worker=None)
    exec(compile(ast.Module(body=selected, type_ignores=[]), str(SOURCE), "exec"), namespace)
    return {name: namespace[name] for name in FUNCTIONS}


def task_probe(name, labels, *, seed=None):
    go = name == FUNCTIONS[0]
    scores = {label: [0.90, 0.95, 1.0] for label in labels}
    scores["A:low"] = [0.01, 0.02, 0.03]
    entries = {label.split(":", 1)[-1]: {
        "interpro_entries": [("IPR", "Repeat", None, "fixture",
            [("GO:1", "molecular_function", "fixture")] if go else [])]}
        for label in [*labels, "A:low"]}
    normalizer = SimpleNamespace(query_interpro_entry_vep_anno=lambda key, _mapping: entries.get(key))
    supplied = set([*labels, "A:low", "A:missing_scores", "A:missing_mapping"])
    scores["A:missing_mapping"] = [0.99]
    kwargs = {} if seed is None else {"seed": seed}
    try:
        procedures()[name](supplied, scores, [0.8, 0.9], normalizer, {}, **kwargs)
    except Submitted as result:
        return [[task[0], task[3], task[4]] for task in result.tasks]
    raise AssertionError("Expected eligible candidate tasks")


def hash_probe():
    labels = [f"G{i}:Pfam:PF{i}" for i in range(200)]
    return {name: task_probe(name, labels) for name in FUNCTIONS}


class RescueOrderTests(unittest.TestCase):
    def test_original_installer_source_stays_unsorted(self):
        before = SOURCE.read_bytes()
        tree = ast.parse(before)
        for name in FUNCTIONS:
            function = next(n for n in tree.body if isinstance(n, ast.FunctionDef) and n.name == name)
            loops = [n for n in ast.walk(function) if isinstance(n, ast.For)
                     and isinstance(n.target, ast.Name) and n.target.id == 'domain_path'
                     and isinstance(n.iter, ast.Name) and n.iter.id == 'tolerant_domains']
            self.assertEqual(len(loops), 1)
        rendered = render_source(before)
        self.assertEqual(SOURCE.read_bytes(), before)
        self.assertEqual(rendered.replace(b'sorted(tolerant_domains)', b'tolerant_domains'), before)

    def test_hash_seed_does_not_change_domain_to_seed(self):
        outcomes = []
        for seed in ("0", "1", "42"):
            result = subprocess.run([sys.executable, str(Path(__file__).resolve()), "--probe"],
                                    env=dict(os.environ, PYTHONHASHSEED=seed),
                                    check=True, capture_output=True, text=True)
            outcomes.append(json.loads(result.stdout))
        self.assertEqual(outcomes[0], outcomes[1])
        self.assertEqual(outcomes[0], outcomes[2])
        for name, initial in zip(FUNCTIONS, (42, 99)):
            self.assertEqual([row[0] for row in outcomes[0][name]],
                             sorted(f"G{i}:Pfam:PF{i}" for i in range(200)))
            self.assertEqual([row[1:] for row in outcomes[0][name]],
                             [[500, initial + i] for i in range(200)])

    def test_unicode_order_and_wrapped_seeds(self):
        labels = ["\U0001f600:domain", "\u00e9:domain", "a:domain", "Z:domain"]
        for name in FUNCTIONS:
            forward = task_probe(name, labels, seed=0xFFFFFFFE)
            reverse = task_probe(name, list(reversed(labels)), seed=0xFFFFFFFE)
            self.assertEqual(forward, reverse)
            self.assertEqual([row[0] for row in forward],
                             ["Z:domain", "a:domain", "\u00e9:domain", "\U0001f600:domain"])
            self.assertEqual([row[2] for row in forward], [0xFFFFFFFE, 0xFFFFFFFF, 0, 1])

    def test_empty_candidates_keep_original_return(self):
        normalizer = SimpleNamespace(query_interpro_entry_vep_anno=lambda *_: None)
        for function in procedures().values():
            self.assertEqual(function(set(), {}, [], normalizer, {}), (set(), {}))


if __name__ == "__main__":
    if sys.argv[1:] == ["--probe"]:
        print(json.dumps(hash_probe()))
    else:
        unittest.main()
