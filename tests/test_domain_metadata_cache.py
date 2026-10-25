"""Real tiny-VCF metadata extraction, cache association and publication faults.

Only actual extractor functions/cache block run; no DAS model, installer, large
input, network or production cache is used. I/O faults are explicitly injected.
The frozen former block establishes successful default-protocol pickle bytes.
"""
import ast
import importlib.util
import json
import logging
import multiprocessing as mp
import os
from pathlib import Path
import pickle
import stat
import subprocess
import sys
import tempfile
import unittest
from unittest import mock

import pysam


ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT / "scripts/am_pick_intolerant_domains.py"
SPEC = importlib.util.spec_from_file_location("domain_metadata_cache_subject", ROOT / "scripts/domain_metadata_cache.py")
CACHE = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(CACHE)
TREE = ast.parse(SOURCE.read_text(), filename=str(SOURCE))
FUNCTIONS = {node.name: node for node in TREE.body if isinstance(node, ast.FunctionDef)}
EXTRACTORS = [FUNCTIONS[name] for name in ("extract_domain_metadata", "_extract_metadata_worker",
                                         "_merge_domain_metadata", "extract_domain_metadata_parallel")]
STAT_SOURCE = ROOT / "scripts/stat_protein_domain_amscores.py"
STAT_TREE = ast.parse(STAT_SOURCE.read_text(), filename=str(STAT_SOURCE))
EXCLUSIONS = [node for node in STAT_TREE.body if isinstance(node, ast.Assign)
              and any(isinstance(item, ast.Name) and item.id == "_NON_DOMAIN_SOURCES" for item in node.targets)]


def compile_nodes(nodes, filename):
    return compile(ast.fix_missing_locations(ast.Module(body=nodes, type_ignores=[])), filename, "exec")


class DomainMetadataCache(unittest.TestCase):
    def setUp(self):
        retained = os.environ.get("PRIVA_DOMAIN_METADATA_FIXTURE_ROOT")
        if retained:
            self.folder = Path(retained) / self._testMethodName
            self.folder.mkdir(parents=True, exist_ok=False)
        else:
            temporary = tempfile.TemporaryDirectory(prefix="domain_metadata_fixture_", dir=os.environ.get(
                "PRIVA_TEST_TMPDIR", os.environ.get("TMPDIR")))
            self.addCleanup(temporary.cleanup)
            self.folder = Path(temporary.name)
        self.cache = self.folder / "domain_metadata.hg38.pkl"
        self.a = self.folder / "synthetic_A.vcf"
        self.b = self.folder / "synthetic_B.vcf"
        self.write_vcf(self.a, "A")
        self.write_vcf(self.b, "B")
        self.log = logging.getLogger("domain_metadata_fixture")
        self.log.addHandler(logging.NullHandler())
        self.namespace = {"os": os, "pickle": pickle, "pysam": pysam, "mp": mp,
                          "logger": self.log, "__file__": str(SOURCE),
                          "self_directory": str(SOURCE.parent),
                          "load_or_extract_domain_metadata": CACHE.load_or_extract_domain_metadata}
        exec(compile_nodes(EXCLUSIONS, str(STAT_SOURCE)), self.namespace)
        exec(compile_nodes(EXTRACTORS, str(SOURCE)), self.namespace)
        self.actual_extract = self.namespace["extract_domain_metadata_parallel"]
        self.calls = []
        self.sources = (SOURCE, STAT_SOURCE)

    def write_vcf(self, path, profile):
        positions, length = ((1, 400), 400) if profile == "A" else ((100, 101), 800)
        text = ('##fileformat=VCFv4.2\n##contig=<ID=chr1,length=1000>\n'
                '##INFO=<ID=CSQ,Number=.,Type=String,Description="Synthetic fixture. Format: Allele|DOMAINS|Gene|SYMBOL|Protein_position|Feature_type|Feature|Unused">\n'
                '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n')
        for number, aa in enumerate(positions, start=1):
            text += (f'chr1\t{number}\t.\tA\tC\t.\t.\t'
                     f'CSQ=C|Pfam:PF00001|ENSG_SYNTHETIC|SYNTHETIC|{aa}/{length}|Transcript|ENST_SYNTHETIC|x\n')
        path.write_text(text)

    def extract(self, source):
        self.calls.append(str(source))
        return self.actual_extract(str(source), threads=1)

    def run_cache(self, source=None, *, extractor=None, assembly="hg38", sources=None):
        selected = source if source is not None else self.a
        return CACHE.load_or_extract_domain_metadata(cache_path=self.cache, vcf_path=selected,
            assembly=assembly, extractor=extractor or (lambda: self.extract(selected)),
            source_paths=sources if sources is not None else self.sources, logger=self.log)

    def snapshot(self, path):
        if not path.exists():
            return None
        value = path.stat()
        return (path.read_bytes(), value.st_ino, value.st_size, value.st_mtime_ns,
                value.st_ctime_ns, value.st_mode, value.st_nlink)

    def pair(self):
        return self.snapshot(self.cache), self.snapshot(CACHE.receipt_path(self.cache))

    def expected(self, source=None, *, assembly="hg38", sources=None):
        return CACHE.input_state(vcf_path=source or self.a, assembly=assembly,
                                 source_paths=sources if sources is not None else self.sources)

    def ready(self, source=None, **kwargs):
        return CACHE.association_matches(cache_path=self.cache, expected_inputs=self.expected(source, **kwargs))

    def fault_case(self, label, prior):
        self.cache = self.folder / label / "domain_metadata.hg38.pkl"
        self.cache.parent.mkdir()
        if prior:
            self.run_cache(self.a)
        source = self.b if prior else self.a
        before = self.pair()
        marker = self.cache.parent / ("." + self.cache.name + ".write-retained")
        marker.mkdir()
        (marker / "older-attempt").write_text("Unrelated older diagnostic fixture")
        return source, before, marker

    def assert_own_staging_cleaned(self, marker):
        self.assertEqual(list(self.cache.parent.glob("." + self.cache.name + ".write-*")), [marker])
        self.assertEqual((marker / "older-attempt").read_text(), "Unrelated older diagnostic fixture")

    def test_real_block_success_bytes_order_and_protocol_match_frozen_before(self):
        legacy_output = self.folder / "legacy"
        legacy_output.mkdir()
        namespace = dict(self.namespace, output_dir=str(legacy_output), assembly="hg38",
                         am_vcf=str(self.a), threads=1)
        baseline = ROOT / "tests/fixtures/das_metadata_cache_before_association_20261018.py"
        exec(compile(baseline.read_text(), str(baseline), "exec"), namespace)
        analyze = FUNCTIONS["analyze_domain_data"]
        start = next(index for index, node in enumerate(analyze.body) if isinstance(node, ast.Assign)
                     and any(isinstance(item, ast.Name) and item.id == "domain_metadata" for item in node.targets))
        namespace = dict(self.namespace, output_dir=str(self.folder), assembly="hg38",
                         am_vcf=str(self.a), threads=1)
        exec(compile_nodes(analyze.body[start:start + 2], str(SOURCE)), namespace)
        expected = (legacy_output / self.cache.name).read_bytes()
        self.assertEqual(self.cache.read_bytes(), expected)
        self.assertEqual(self.cache.read_bytes()[:2], bytes((0x80, pickle.DEFAULT_PROTOCOL)))
        with self.cache.open("rb") as handle:
            loaded = pickle.load(handle)
        self.assertEqual(list(loaded), ["ENSG_SYNTHETIC:Pfam:PF00001"])
        self.assertEqual(list(loaded.values())[0], {
            "gene_symbol": "SYNTHETIC", "aa_start": 1, "aa_end": 400, "prot_length": 400})
        self.assertTrue(self.ready())

    def test_matching_receipt_warm_reuse_preserves_bytes_and_identity(self):
        expected = self.run_cache()
        before = self.pair()
        self.assertEqual(self.run_cache(extractor=mock.Mock(side_effect=AssertionError("unexpected extraction"))), expected)
        self.assertEqual(self.pair(), before)

    def test_selected_vcf_change_during_warm_load_fails_without_cache_mutation(self):
        self.run_cache()
        before = self.pair()
        actual_load = pickle.load

        def changing_input(handle):
            value = actual_load(handle)
            self.write_vcf(self.a, "B")
            return value

        with mock.patch.object(pickle, "load", side_effect=changing_input):
            with self.assertRaisesRegex(ValueError, "changed during cache loading"):
                self.run_cache()
        self.assertEqual(self.pair(), before)
        self.assertEqual(len(self.calls), 1)

    def test_output_replacement_during_warm_load_fails_without_rewriting(self):
        self.run_cache()
        old_receipt = self.snapshot(CACHE.receipt_path(self.cache))
        actual_load = pickle.load
        replacement_snapshot = []

        def changing_output(handle):
            value = actual_load(handle)
            replacement = self.folder / "externally_replaced_cache.pkl"
            replacement.write_bytes(pickle.dumps(value))
            os.replace(replacement, self.cache)
            replacement_snapshot.append(self.snapshot(self.cache))
            return value

        with mock.patch.object(pickle, "load", side_effect=changing_output):
            with self.assertRaisesRegex(ValueError, "changed during cache loading"):
                self.run_cache()
        self.assertEqual(self.snapshot(self.cache), replacement_snapshot[0])
        self.assertEqual(self.snapshot(CACHE.receipt_path(self.cache)), old_receipt)
        self.assertEqual(len(self.calls), 1)

    def test_unreceipted_valid_old_profile_regenerates_selected_input(self):
        self.cache.write_bytes(pickle.dumps(self.actual_extract(str(self.a), threads=1)))
        old = self.snapshot(self.cache)

        def selected():
            self.assertEqual(self.snapshot(self.cache), old)
            return self.extract(self.b)

        result = self.run_cache(self.b, extractor=selected)
        self.assertEqual(result, self.actual_extract(str(self.b), threads=1))
        self.assertEqual(self.calls, [str(self.b)])
        self.assertTrue(self.ready(self.b))

    def test_unreceipted_partial_pickle_recovers_from_selected_vcf(self):
        self.cache.write_bytes(b"\x80\x04\x95\x20")
        self.assertEqual(self.run_cache(), self.actual_extract(str(self.a), threads=1))
        self.assertTrue(self.ready())

    def test_different_selected_path_regenerates_current_profile(self):
        first = self.run_cache(self.a)
        result = self.run_cache(self.b)
        self.assertNotEqual(first, result)
        self.assertEqual(result, self.actual_extract(str(self.b), threads=1))
        self.assertEqual(self.calls, [str(self.a), str(self.b)])

    def test_changed_vcf_bytes_and_mtime_at_same_path_regenerate(self):
        first = self.run_cache()
        before = self.a.stat()
        self.write_vcf(self.a, "B")
        os.utime(self.a, ns=(before.st_atime_ns, before.st_mtime_ns + 10_000_000_000))
        result = self.run_cache()
        self.assertNotEqual(first, result)
        self.assertEqual(result, self.actual_extract(str(self.a), threads=1))
        self.assertEqual(self.calls, [str(self.a), str(self.a)])

    def test_assembly_alias_reuses_and_other_assembly_invalidates(self):
        self.run_cache(assembly="hg38")
        before = self.pair()
        self.run_cache(assembly="GRCh38")
        self.assertEqual(self.pair(), before)
        self.assertEqual(len(self.calls), 1)
        self.run_cache(assembly="hg19")
        self.assertEqual(len(self.calls), 2)
        self.assertTrue(self.ready(assembly="GRCh37"))

    def test_extraction_source_change_invalidates_association(self):
        implementation = self.folder / "synthetic_extractor_implementation.py"
        implementation.write_text("# version one\n")
        sources = (*self.sources, implementation)
        self.run_cache(sources=sources)
        implementation.write_text("# version two changed\n")
        self.run_cache(sources=sources)
        self.assertEqual(len(self.calls), 2)
        self.assertTrue(self.ready(sources=sources))

    def test_vcf_change_during_generation_preserves_prior_pair(self):
        self.run_cache()
        before = self.pair()

        def changing_input():
            result = self.extract(self.b)
            self.write_vcf(self.b, "A")
            return result

        with self.assertRaisesRegex(ValueError, "changed during"):
            self.run_cache(self.b, extractor=changing_input)
        self.assertEqual(self.pair(), before)

    def test_source_change_during_generation_preserves_prior_pair(self):
        implementation = self.folder / "synthetic_extractor_implementation.py"
        implementation.write_text("# version one\n")
        sources = (*self.sources, implementation)
        self.run_cache(sources=sources)
        before = self.pair()

        def changing_source():
            result = self.extract(self.b)
            implementation.write_text("# version changed during extraction\n")
            return result

        with self.assertRaisesRegex(ValueError, "changed during"):
            self.run_cache(self.b, extractor=changing_source, sources=sources)
        self.assertEqual(self.pair(), before)

    def test_partial_dump_failure_preserves_prior_and_retry_completes(self):
        for prior in (False, True):
            with self.subTest(prior=prior):
                source, before, marker = self.fault_case("dump-" + str(prior), prior)

                def partial_dump(value, handle):
                    handle.write(b"\x80\x04\x95\x20")
                    handle.flush()
                    raise OSError("Injected partial pickle dump")

                with mock.patch.object(pickle, "dump", side_effect=partial_dump):
                    with self.assertRaises(OSError):
                        self.run_cache(source)
                self.assertEqual(self.pair(), before)
                self.assert_own_staging_cleaned(marker)
                self.run_cache(source)
                self.assertTrue(self.ready(source))

    def write_io_failure(self, operation):
        for prior in (False, True):
            with self.subTest(prior=prior, operation=operation):
                source, before, marker = self.fault_case(operation + "-" + str(prior), prior)
                actual_open = Path.open

                class FileBoundary:
                    def __init__(self, original):
                        self.original = original

                    def __getattr__(self, name):
                        return getattr(self.original, name)

                    def __enter__(self):
                        self.original.__enter__()
                        return self

                    def flush(self):
                        if operation == "flush":
                            raise OSError("Injected flush failure")
                        return self.original.flush()

                    def __exit__(self, *args):
                        result = self.original.__exit__(*args)
                        if operation == "close":
                            raise OSError("Injected close failure")
                        return result

                def staged_open(path, *args, **kwargs):
                    result = actual_open(path, *args, **kwargs)
                    if path.name == self.cache.name and path.parent.name.startswith("." + self.cache.name + ".write-"):
                        return FileBoundary(result)
                    return result

                with mock.patch.object(Path, "open", staged_open):
                    with self.assertRaises(OSError):
                        self.run_cache(source)
                self.assertEqual(self.pair(), before)
                self.assert_own_staging_cleaned(marker)

    def test_flush_failure_preserves_prior_pair(self):
        self.write_io_failure("flush")

    def test_close_failure_preserves_prior_pair(self):
        self.write_io_failure("close")

    def test_data_replace_failure_preserves_prior_pair(self):
        for prior in (False, True):
            with self.subTest(prior=prior):
                source, before, marker = self.fault_case("data-replace-" + str(prior), prior)
                replace = os.replace

                def failing_replace(first, second):
                    if Path(second) == self.cache:
                        raise OSError("Injected data replace failure")
                    return replace(first, second)

                with mock.patch.object(os, "replace", side_effect=failing_replace):
                    with self.assertRaises(OSError):
                        self.run_cache(source)
                self.assertEqual(self.pair(), before)
                self.assert_own_staging_cleaned(marker)

    def receipt_failure(self, operation):
        for prior in (False, True):
            with self.subTest(prior=prior):
                source, before, marker = self.fault_case("receipt-" + operation + "-" + str(prior), prior)
                sidecar = CACHE.receipt_path(self.cache)
                replace = os.replace

                def failing_replace(first, second):
                    if Path(second) == sidecar:
                        raise OSError("Injected receipt replace failure")
                    return replace(first, second)

                def failing_json(value, handle, **kwargs):
                    handle.write('{"version":')
                    raise OSError("Injected receipt write failure")

                boundary = mock.patch.object(os, "replace", side_effect=failing_replace) if operation == "replace" else \
                    mock.patch.object(json, "dump", side_effect=failing_json)
                with boundary:
                    with self.assertRaises(OSError):
                        self.run_cache(source)
                self.assertEqual(self.snapshot(sidecar), before[1])
                self.assertFalse(self.ready(source))
                self.assert_own_staging_cleaned(marker)
                calls = len(self.calls)
                self.run_cache(source)
                self.assertEqual(len(self.calls), calls + 1)
                self.assertTrue(self.ready(source))

    def test_receipt_write_failure_is_not_ready_and_retry_completes(self):
        self.receipt_failure("write")

    def test_receipt_replace_failure_is_not_ready_and_retry_completes(self):
        self.receipt_failure("replace")

    def test_changed_output_metadata_regenerates_without_loading_partial(self):
        self.run_cache()
        self.cache.write_bytes(b"\x80\x04\x95\x20")
        self.run_cache()
        self.assertEqual(len(self.calls), 2)
        self.assertTrue(self.ready())

    def test_original_missing_prot_length_rule_still_regenerates(self):
        self.run_cache(extractor=lambda: {"legacy": {"aa_start": 1, "aa_end": 2}})
        self.assertTrue(self.ready())
        result = self.run_cache()
        self.assertIn("ENSG_SYNTHETIC:Pfam:PF00001", result)
        self.assertEqual(len(self.calls), 1)

    def test_identity_is_device_independent_and_does_not_hash_vcf(self):
        original = Path.read_bytes

        def source_only(path):
            if path == self.a:
                raise AssertionError("VCF bytes must not be hashed for association")
            return original(path)

        with mock.patch.object(Path, "read_bytes", source_only):
            value = self.expected()
        self.assertEqual(set(value["vcf"]), {"path", "target", "inode", "size", "mtime_ns"})
        for description in [value["vcf"], *value["implementation"]]:
            self.assertNotIn("device", description)
            self.assertNotIn("st_dev", description)

    def test_existing_mode_and_new_file_umask_are_preserved(self):
        self.run_cache()
        os.chmod(self.cache, 0o640)
        self.run_cache(self.b)
        self.assertEqual(stat.S_IMODE(self.cache.stat().st_mode), 0o640)
        self.cache = self.folder / "new-cache.pkl"
        old = os.umask(0o027)
        try:
            self.run_cache()
        finally:
            os.umask(old)
        self.assertEqual(stat.S_IMODE(self.cache.stat().st_mode), 0o640)

    def test_missing_selected_vcf_fails_without_replacement_or_prior_mutation(self):
        self.run_cache()
        before = self.pair()
        with self.assertRaises(FileNotFoundError):
            self.run_cache(self.folder / "missing-selected.vcf")
        self.assertEqual(self.pair(), before)

    def test_named_cli_checks_metadata_only_and_rejects_mismatch(self):
        self.run_cache()
        before = self.pair()
        command = [sys.executable, str(ROOT / "scripts/domain_metadata_cache.py"), "check",
                   "--cache", str(self.cache), "--vcf", str(self.a), "--assembly", "GRCh38"]
        for source in self.sources:
            command += ["--source-file", str(source)]
        result = subprocess.run(command, capture_output=True, text=True, timeout=10)
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertEqual(result.stdout, "")
        command[command.index(str(self.a))] = str(self.b)
        result = subprocess.run(command, capture_output=True, text=True, timeout=10)
        self.assertNotEqual(result.returncode, 0)
        self.assertEqual(self.pair(), before)


if __name__ == "__main__":
    unittest.main()
