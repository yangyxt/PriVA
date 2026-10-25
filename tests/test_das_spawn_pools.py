"""Real small-model fork/spawn comparisons for the five independent DAS pools.

Run only in a bounded compute allocation; keep generated payloads outside source:
  python test_das_spawn_pools.py prepare --repository-root REPO \
      --baseline-source FROZEN_57cc_SOURCE --output-root FIXTURE_ROOT
  python FIXTURE_ROOT/source/test_das_spawn_pools.py run \
      --work-root FIXTURE_ROOT --pins-sha256 PREPARE_PRINTED_SHA256

The baseline is the actual pre-fix source, not an independently rewritten oracle.
Only the main KS orchestration is extracted by AST; its worker/statistics and all
other tested entry points execute unchanged. Completion orders and files remain
recorded as produced. Per-domain comparisons do not assert that unordered result
tables, set iteration or a full scientific run are deterministic.
"""
import argparse
import ast
from contextlib import contextmanager
import copy
import hashlib
import importlib
import importlib.util
import io
import json
import multiprocessing as real_mp
import os
from pathlib import Path
import pickle
import resource
import shutil
import sys
import unittest


BASELINE_SHA256 = "57cc7be79168bcbc10aebc12efe94ea019fd33c148f35b5ba65360fc29b648c3"
CONTEXT_COMMENT = "    # Spawn independent tasks without inheriting the genome-wide score model.\n"
DEPENDENCIES = ("am_pick_intolerant_domains.py", "protein_domain_mapping.py",
                "interpro_mapping_io.py", "domain_metadata_cache.py",
                "stat_protein_domain_amscores.py")
CONFIG = None


def sha256(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def freeze(value):
    """Exact typed evidence, retaining dict/list order and binary float values."""
    import numpy as np
    if isinstance(value, np.ndarray):
        return ["ndarray", value.dtype.str, list(value.shape), value.tobytes().hex(),
                bool(value.flags.writeable)]
    if isinstance(value, np.generic):
        return ["numpy", value.dtype.str, value.tobytes().hex()]
    if isinstance(value, float):
        return ["float", value.hex()]
    if isinstance(value, dict):
        return ["dict", [[freeze(k), freeze(v)] for k, v in value.items()]]
    if isinstance(value, (tuple, list)):
        return [type(value).__name__, [freeze(v) for v in value]]
    if isinstance(value, set):
        # A comparison representation only; never fed back to scientific code.
        return ["set", sorted(value)]
    if value is None or isinstance(value, (str, int, bool)):
        return value
    raise TypeError(f"Unsupported fixture evidence type: {type(value).__name__}")


def runtime_probe(module_name):
    module = importlib.import_module(module_name)
    return {"start_method": real_mp.get_start_method(), "pid": os.getpid(),
            "inherited_parent_model": hasattr(module, "_fixture_parent_model"),
            "max_rss_kib": resource.getrusage(resource.RUSAGE_SELF).ru_maxrss}


class TracedPool:
    """Observe actual pools without replacing workers, payloads or iterator APIs."""
    def __init__(self, owner, method, processes):
        self.record = {"context": method, "processes": processes}
        owner.records.append(self.record)
        self.pool = real_mp.get_context(method).Pool(processes)
        self.record["probe"] = self.pool.apply(runtime_probe, (owner.module.__name__,))
        self.module = owner.module

    def __enter__(self):
        self.pool.__enter__()
        return self

    def __exit__(self, exc_type, exc_value, traceback):
        self.record["exit_exception"] = None if exc_type is None else exc_type.__name__
        return self.pool.__exit__(exc_type, exc_value, traceback)

    def observe(self, worker, tasks, api):
        assert isinstance(tasks, list), "The original call sites construct lists"
        self.worker = worker
        self.tasks = copy.deepcopy(tasks)
        self.outputs = []
        self.record.update(api=api, worker=worker.__name__, completion_order=[])
        semantic_tasks = [task[:2] if worker.__name__ == "visualize_single_domain" else task
                          for task in tasks]
        self.record["payloads"] = freeze(semantic_tasks)
        if worker.__name__ in ("process_domain_ks", "_dip_rescue_worker"):
            self.record["seeds"] = [task[-1] for task in tasks]

    def remember(self, result):
        self.outputs.append(copy.deepcopy(result))
        if isinstance(result, tuple):
            key = result[0]
        elif isinstance(result, dict) and "domain" in result:
            key = result["domain"]
        else:
            key = None
        self.record["completion_order"].append(key)
        self.record["returns"] = freeze(self.outputs)

    def map(self, worker, tasks):
        self.observe(worker, tasks, "map")
        results = self.pool.map(worker, tasks)
        for result in results:
            self.remember(result)
        return results

    def imap_unordered(self, worker, tasks):
        self.observe(worker, tasks, "imap_unordered")
        for result in self.pool.imap_unordered(worker, tasks):
            self.remember(result)
            yield result


class TracedMP:
    def __init__(self, module):
        self.module, self.records, self.pools = module, [], []

    def Pool(self, processes):
        pool = TracedPool(self, "fork", processes)
        self.pools.append(pool)
        return pool

    def get_context(self, method):
        assert method == "spawn", method
        owner = self

        class Context:
            def Pool(self, processes):
                pool = TracedPool(owner, method, processes)
                owner.pools.append(pool)
                return pool
        return Context()


@contextmanager
def traced(module):
    original = module.mp
    trace = TracedMP(module)
    module.mp = trace
    try:
        yield trace
    finally:
        module.mp = original


class FixtureNormalizer:
    def __init__(self, entries):
        self.entries = entries

    def query_interpro_entry_vep_anno(self, name, mapping):
        return self.entries.get(name)


def entry(kind, go_mf):
    return {"interpro_entries": [("IPR_FIXTURE", kind, None, "fixture",
                                   [("GO:fixture", "molecular_function", "fixture")] if go_mf else [], "fixture")]}


def load_module(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


def ks_block(module, output_dir, *, cap=60):
    """Execute the actual main KS task construction, pool and result-table block."""
    import numpy as np
    tree = ast.parse(Path(module.__file__).read_text())
    analyze = next(node for node in tree.body if isinstance(node, ast.FunctionDef)
                   and node.name == "analyze_domain_data")

    def assignment(node, name):
        return isinstance(node, ast.Assign) and any(isinstance(t, ast.Name) and t.id == name for t in node.targets)

    start = next(i for i, node in enumerate(analyze.body) if assignment(node, "cap"))
    end = next(i for i, node in enumerate(analyze.body) if assignment(node, "not_more_tolerant"))
    wrapper = ast.parse("def _fixture_main_ks():\n    pass\n").body[0]
    wrapper.body = analyze.body[start:end] + [ast.Return(value=ast.Call(func=ast.Name(id="locals", ctx=ast.Load()),
                                                                     args=[], keywords=[]))]
    block = compile(ast.fix_missing_locations(ast.Module(body=[wrapper], type_ignores=[])), module.__file__, "exec")
    inputs = {"query_domains": [("query:fallback", np.linspace(.05, .4, 9)),
                                ("query:typed", np.linspace(.4, .9, 12)),
                                ("query:short", np.array([.8, .95]))],
              "query_type_map": {"query:typed": "Family"},
              "type_composites": {"Family": np.linspace(.2, .99, 20)},
              "fixed_composite": np.linspace(.1, .99, 75)}
    before = freeze(inputs)
    namespace = dict(module.__dict__, **inputs, T=cap, threads=2, seed_global=42,
                     output_dir=str(output_dir), assembly="hg38", fdr_threshold=.01)
    exec(block, namespace)
    namespace.update(namespace["_fixture_main_ks"]())
    assert freeze(inputs) == before, "Main KS inputs mutated"
    return namespace["ks_results"], namespace["results_df"]


class DASPoolComparisons(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        if CONFIG is None:
            raise unittest.SkipTest("Use the named prepare/run commands in a bounded Slurm allocation")
        cls.root = CONFIG
        sys.path[:0] = [str(cls.root / "source/scripts"), str(cls.root / "source")]
        cls.baseline = load_module("das_fork_baseline", cls.root / "source/das_fork_baseline.py")
        cls.candidate = load_module("am_pick_intolerant_domains", cls.root / "source/scripts/am_pick_intolerant_domains.py")
        # This stands for a loaded parent-only model, not a full-cache RAM claim.
        cls.baseline._fixture_parent_model = bytearray(8 * 1024 * 1024)
        cls.candidate._fixture_parent_model = bytearray(8 * 1024 * 1024)

    def setUp(self):
        self.folder = self.root / "fixtures" / self._testMethodName
        self.folder.mkdir(parents=True, exist_ok=False)
        self.traces = []

    def tearDown(self):
        (self.folder / "pool-traces.json").write_text(json.dumps(self.traces, indent=2) + "\n")

    def pair(self, operation, *, expected_api="imap_unordered", failure=None):
        results, traces = [], []
        for name, module in (("fork", self.baseline), ("spawn", self.candidate)):
            folder = self.folder / name
            folder.mkdir()
            with traced(module) as trace:
                if failure:
                    with self.assertRaises(failure) as error:
                        operation(module, folder)
                    results.append((type(error.exception).__name__, str(error.exception)))
                else:
                    results.append(operation(module, folder))
            self.traces.append({"variant": name, "pools": trace.records})
            traces.append(trace)
            self.assertEqual(len(trace.records), 1)
            record = trace.records[0]
            self.assertEqual(record["context"], name)
            self.assertEqual(record["probe"]["start_method"], name)
            self.assertEqual(record["probe"]["inherited_parent_model"], name == "fork")
            self.assertEqual(record["api"], expected_api)
            self.assertEqual(record["exit_exception"] is not None, bool(failure))
        self.assertEqual(traces[0].records[0]["payloads"], traces[1].records[0]["payloads"])
        if failure:
            self.assertEqual(results[0], results[1])
            return results, traces
        first, second = (trace.pools[0] for trace in traces)
        if expected_api == "map":
            self.assertEqual(freeze(first.outputs), freeze(second.outputs))
        else:
            def keyed(outputs):
                return {r[0] if isinstance(r, tuple) else r["domain"]: freeze(r)
                        for r in outputs if r is not None}
            self.assertEqual(keyed(first.outputs), keyed(second.outputs))
            self.assertEqual(sum(r is None for r in first.outputs), sum(r is None for r in second.outputs))
        for pool in (first, second):
            if pool.worker.__name__ != "visualize_single_domain":
                serial = [pool.worker(task) for task in pool.tasks]
                actual = pool.outputs
                if expected_api == "map":
                    self.assertEqual(freeze(serial), freeze(actual))
                else:
                    self.assertEqual({r[0]: freeze(r) for r in serial}, {r[0]: freeze(r) for r in actual})
        return results, traces

    def test_only_five_contexts_and_explanatory_comment_changed(self):
        baseline = Path(self.baseline.__file__).read_text()
        candidate = Path(self.candidate.__file__).read_text()
        self.assertEqual(hashlib.sha256(baseline.encode()).hexdigest(), BASELINE_SHA256)
        self.assertEqual(candidate.count("mp.get_context('spawn').Pool("), 5)
        self.assertEqual(candidate.replace("mp.get_context('spawn').Pool(", "mp.Pool(")
                         .replace(CONTEXT_COMMENT, ""), baseline)
        tree = ast.parse(candidate)
        owners = []
        for node in tree.body:
            if isinstance(node, ast.FunctionDef):
                owners.extend([node.name for item in ast.walk(node) if isinstance(item, ast.Call)
                               and ast.unparse(item.func) == "mp.get_context('spawn').Pool"])
        self.assertEqual(owners, ["visualize_domain_distribution", "rescue_bimodal_go_mf_domains",
                                 "rescue_non_gomf_domains", "extract_domain_metadata_parallel", "analyze_domain_data"])

    def test_visualization_real_png_statistics_skips_errors_and_nonmutation(self):
        import numpy as np
        scores = {"ENSG_FIXTURE": {"Pfam": {
            "A": {"max_distribution": np.linspace(.05, .95, 8)},
            "B": {"max_distribution": np.array([.05, .08, .2, .3, .5, .6, .8, .81, .9, .95, .99])},
            "SIX": {"max_distribution": np.linspace(.1, .8, 6)},
            "SEVEN": {"max_distribution": np.linspace(.1, .8, 7)},
            "EXISTS": {"max_distribution": np.linspace(.1, .8, 8)},
            "BAD": {"max_distribution": np.array(["invalid"] * 8)},
            "BELOW": {"max_distribution": np.linspace(.1, .8, 5)}}}}
        before = freeze(scores)

        def operation(module, folder):
            target = folder / "am_domain_distributions"
            target.mkdir()
            existing = target / "distribution_ENSG_FIXTURE:Pfam:EXISTS.png"
            existing.write_bytes(b"preserved fixture sentinel")
            result = module.visualize_domain_distribution(scores, str(folder), threads=2, assembly="hg38")
            self.assertEqual(freeze(scores), before)
            self.assertEqual(existing.read_bytes(), b"preserved fixture sentinel")
            self.assertEqual(list(result), [f"ENSG_FIXTURE:Pfam:{key}" for key in ("A", "B", "SIX", "SEVEN", "EXISTS", "BAD")])
            return result

        results, _ = self.pair(operation)
        self.assertEqual(freeze(results[0]), freeze(results[1]))
        for suffix in ("distribution_ENSG_FIXTURE:Pfam:A.png", "distribution_ENSG_FIXTURE:Pfam:B.png", "domain_amscores.hg38.pkl"):
            files = [self.folder / name / "am_domain_distributions" / suffix for name in ("fork", "spawn")]
            self.assertEqual(files[0].read_bytes(), files[1].read_bytes())
            self.assertGreater(files[0].stat().st_size, 0)
        self.compare_tables("am_domain_distributions/distribution_test_results.hg38.tsv", "domain")

    def compare_tables(self, relative, key):
        import pandas as pd
        tables = [pd.read_table(self.folder / name / relative) for name in ("fork", "spawn")]
        self.assertEqual(list(tables[0].columns), list(tables[1].columns))
        # Keep native table bytes/order. Compare corresponding task rows only.
        self.assertEqual({row[key]: freeze(row) for row in tables[0].to_dict("records")},
                         {row[key]: freeze(row) for row in tables[1].to_dict("records")})

    def rescue_case(self, go_mf, *, seed=None, fail=False):
        import numpy as np
        prefix = "GO" if go_mf else "NON"
        paths = [f"ENSG_FIXTURE:Pfam:{prefix}_{i}" for i in range(3)]
        scores = {paths[0]: np.linspace(.65, .95, 9), paths[1]: np.array([.78] * 5 + [.99] * 5),
                  paths[2]: np.linspace(.001, .2, 9)}
        skipped = "ENSG_FIXTURE:Pfam:SKIP"
        scores[skipped] = np.linspace(.8, .99, 8)
        tolerant = set(paths + [skipped, "ENSG_FIXTURE:Pfam:ABSENT"])
        mapping = {p.split(":", 1)[1]: entry("Family" if go_mf else "Repeat", go_mf) for p in paths}
        mapping[skipped.split(":", 1)[1]] = entry("Domain" if go_mf else "Family", go_mf)
        normalizer = FixtureNormalizer(mapping)
        composite = np.linspace(.1, .99, 2 if fail else 80)
        typed = {} if fail else {"Family" if go_mf else "Repeat": np.linspace(.25, .97, 64)}
        inputs_before = freeze((scores, tolerant, composite, typed, mapping))

        def operation(module, folder):
            function = module.rescue_bimodal_go_mf_domains if go_mf else module.rescue_non_gomf_domains
            kwargs = {"seed": seed} if seed is not None else {}
            try:
                return function(tolerant, scores, composite, normalizer, {}, type_composites=typed,
                                threads=2, n_permutations=17, **kwargs)
            finally:
                self.assertEqual(freeze((scores, tolerant, composite, typed, mapping)), inputs_before)

        results, traces = self.pair(operation, failure=ValueError if fail else None)
        if not fail:
            self.assertEqual(freeze(results[0]), freeze(results[1]))
        expected_seed = seed if seed is not None else (42 if go_mf else 99)
        for trace in traces:
            self.assertEqual(trace.records[0]["seeds"], [(expected_seed + i) & 0xFFFFFFFF for i in range(2)])

    def test_go_mf_dip_pool_exact_results_default_seeds_and_nonmutation(self):
        self.rescue_case(True)

    def test_non_go_mf_dip_pool_exact_results_default_seeds_and_nonmutation(self):
        self.rescue_case(False)

    def test_go_mf_existing_seed_wrap_is_preserved(self):
        self.rescue_case(True, seed=0xFFFFFFFF)

    def test_go_mf_dip_worker_failure_propagates(self):
        self.rescue_case(True, fail=True)

    def test_non_go_mf_dip_worker_failure_propagates(self):
        self.rescue_case(False, fail=True)

    def metadata_vcf(self, *, broken=False):
        import pysam
        path = self.folder / "input.vcf"
        fields = "Allele|BROKEN|Gene|SYMBOL|Protein_position|Feature_type|Feature|Unused" if broken else \
                 "Allele|DOMAINS|Gene|SYMBOL|Protein_position|Feature_type|Feature|Unused"
        text = ('##fileformat=VCFv4.2\n##contig=<ID=chr2,length=1000>\n##contig=<ID=chr1,length=1000>\n'
                f'##INFO=<ID=CSQ,Number=.,Type=String,Description="Fixture. Format: {fields}">\n'
                '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n'
                'chr2\t1\t.\tA\tC\t.\t.\tCSQ=C|Pfam:PF01&SMART:SM01|ENSG_FIXTURE|FIRST|20/300|Transcript|ENST_FIXTURE|x\n'
                'chr1\t1\t.\tA\tC\t.\t.\tCSQ=C|Pfam:PF01|ENSG_FIXTURE|LATER|2/400|Transcript|ENST_FIXTURE|x\n')
        path.write_text(text)
        compressed = str(path) + ".gz"
        pysam.tabix_compress(str(path), compressed)
        pysam.tabix_index(compressed, preset="vcf")
        return compressed

    def test_metadata_real_indexed_vcf_order_merge_target_and_nonmutation(self):
        path = self.metadata_vcf()
        target = {"ENSG_FIXTURE:Pfam:PF01"}
        before = (sha256(path), sha256(path + ".tbi"), freeze(target))

        def operation(module, folder):
            result = module.extract_domain_metadata_parallel(path, target_domains=target, threads=2)
            self.assertEqual((sha256(path), sha256(path + ".tbi"), freeze(target)), before)
            return result

        results, _ = self.pair(operation, expected_api="map")
        self.assertEqual(freeze(results[0]), freeze(results[1]))
        self.assertEqual(results[1], {"ENSG_FIXTURE:Pfam:PF01": {
            "gene_symbol": "FIRST", "aa_start": 2, "aa_end": 20, "prot_length": 400}})

    def test_metadata_worker_failure_propagates(self):
        path = self.metadata_vcf(broken=True)
        self.pair(lambda module, folder: module.extract_domain_metadata_parallel(path, threads=2),
                  expected_api="map", failure=ValueError)

    def test_main_ks_actual_task_block_exact_seeds_results_tables_and_nonmutation(self):
        _, traces = self.pair(lambda module, folder: ks_block(module, folder))
        for trace in traces:
            self.assertEqual(trace.records[0]["seeds"], [42, 43, 44])
            self.assertEqual([task[-2] for task in trace.pools[0].tasks], [60, 20, 60])
        self.compare_tables("domain_tolerance_analysis.hg38.tsv", "domain")

    def test_main_ks_worker_failure_propagates(self):
        self.pair(lambda module, folder: ks_block(module, folder, cap=-1), failure=ValueError)


def prepare(*, repository_root, baseline_source, output_root):
    repository_root, baseline_source, output_root = map(Path, (repository_root, baseline_source, output_root))
    if sha256(baseline_source) != BASELINE_SHA256:
        raise ValueError("Expected the exact frozen pre-fix DAS source")
    output_root.mkdir(parents=True, exist_ok=False)
    source = output_root / "source"
    (source / "scripts").mkdir(parents=True)
    files = [(baseline_source, source / "das_fork_baseline.py"),
             (Path(__file__), source / "test_das_spawn_pools.py")]
    files.extend((repository_root / "scripts" / name, source / "scripts" / name) for name in DEPENDENCIES)
    pins = {}
    for original, destination in files:
        before = sha256(original)
        shutil.copyfile(original, destination)
        if sha256(original) != before or sha256(destination) != before:
            raise RuntimeError(f"Source changed during fixture snapshot: {original.name}")
        pins[str(destination.relative_to(output_root))] = before
    manifest = output_root / "source-pins.json"
    manifest.write_text(json.dumps({"schema": "priva.das-spawn-fixture-source.v1", "files": pins}, indent=2) + "\n")
    print(json.dumps({"source_manifest_sha256": sha256(manifest), "files": pins}, sort_keys=True))


def run(*, work_root, pins_sha256):
    global CONFIG
    CONFIG = Path(work_root)
    manifest = CONFIG / "source-pins.json"
    if sha256(manifest) != pins_sha256:
        raise ValueError("Source manifest mismatch")
    pins = json.loads(manifest.read_text())["files"]

    def verify():
        for relative, expected in pins.items():
            if sha256(CONFIG / relative) != expected:
                raise ValueError(f"Frozen source changed: {relative}")

    verify()
    # Explicit fork control is fixture-only; production selects no global method.
    original_default = real_mp.get_context().get_start_method()
    stream = io.StringIO()
    suite = unittest.defaultTestLoader.loadTestsFromTestCase(DASPoolComparisons)
    result = unittest.TextTestRunner(stream=stream, verbosity=2).run(suite)
    (CONFIG / "unittest.log").write_text(stream.getvalue())
    verify()
    if real_mp.get_context().get_start_method() != original_default:
        raise RuntimeError("Global multiprocessing default changed")
    summary = {"schema": "priva.das-spawn-fixture-results.v1", "source_manifest_sha256": pins_sha256,
               "tests": result.testsRun, "failures": len(result.failures), "errors": len(result.errors),
               "skipped": len(result.skipped), "successful": result.wasSuccessful(),
               "global_start_method_before_after": original_default,
               "comparison": "exact per-task typed returns, payloads and seeds; native completion orders retained",
               "scope": "synthetic small models and actual workers; no full-cache acceptance or memory bound",
               "known_original_ordering": "imap_unordered table insertion order and later tolerant-set candidate/seed order are not normalized"}
    (CONFIG / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary, sort_keys=True))
    if not result.wasSuccessful():
        print(stream.getvalue())
        raise SystemExit(1)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    commands = parser.add_subparsers(dest="command", required=True)
    prepare_parser = commands.add_parser("prepare")
    prepare_parser.add_argument("--repository-root", required=True)
    prepare_parser.add_argument("--baseline-source", required=True)
    prepare_parser.add_argument("--output-root", required=True)
    run_parser = commands.add_parser("run")
    run_parser.add_argument("--work-root", required=True)
    run_parser.add_argument("--pins-sha256", required=True)
    args = vars(parser.parse_args())
    command = args.pop("command")
    {"prepare": prepare, "run": run}[command](**args)


if __name__ == "__main__":
    main()
