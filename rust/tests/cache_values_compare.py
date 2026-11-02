#!/usr/bin/env python3
"""Compare every native cache value with CPython's independent pickle decoder.

Run only on trusted existing PriVA pickles: pickle.load can execute code.
  python cache_values_compare.py compare --binary PRIVA --cache NAME \
      --pickle SOURCE --archive CACHE --report REPORT.json --mismatches DIFF.jsonl
  python cache_values_compare.py selftest

The native reader streams tagged values, so its Python decoder is not its own oracle.
Float bits, scalar/container types, duplicates, missingness and meaningful order are
exact. Only the frozen format's fixed-record and lookup-key orders are normalized.
Every unequal leaf/length and missing/extra top key is recorded. No numeric tolerance.
"""
from __future__ import annotations

import argparse
import gzip
import hashlib
import json
import math
from pathlib import Path
import pickle
import resource
import struct
import subprocess
import sys
import time
import unittest


def encode(value: object) -> list:
    """Types stay distinct even where Python's == would equate them."""
    kind = type(value)
    if value is None:
        return ["none"]
    if kind is bool:
        return ["bool", value]
    if kind is int:
        return ["int", str(value)]
    if kind is float:
        return ["float", struct.pack(">d", value).hex()]
    if kind is str:
        return ["str", value]
    if kind is dict:
        return ["dict", [[encode(key), encode(item)] for key, item in value.items()]]
    if kind in (list, tuple, set, frozenset):
        values = [encode(item) for item in value]
        if kind in (set, frozenset):
            values.sort(key=serialized)
        return [kind.__name__, values]
    raise TypeError(f"unsupported reference type: {kind.__module__}.{kind.__name__}")


def serialized(value: object) -> str:
    return json.dumps(value, ensure_ascii=False, separators=(",", ":"), allow_nan=False)


def key_value(tagged: list) -> object:
    if tagged == ["none"]:
        return None
    if tagged[0] == "str" and len(tagged) == 2 and type(tagged[1]) is str:
        return tagged[1]
    if tagged[0] == "tuple" and len(tagged) == 2:
        return tuple(key_value(item) for item in tagged[1])
    raise ValueError("unsupported exported lookup key")


def normalize(cache: str, top_key: object, value: list, path: tuple = ()) -> list:
    """Normalize only ordering explicitly discarded by CACHE_FORMATS.md.

    Record structs do not store their field order. Inner free-key dictionaries do.
    PM1's two transcript lookup maps explicitly sort keys in their frozen schema.
    """
    kind = value[0]
    if kind == "dict":
        pairs = [[key, normalize(cache, top_key, item, path + (key_value(key),))]
                 for key, item in value[1]]
        record = ((cache == "clinvar_aa_stat" and len(path) == 2)
                  or (cache == "clinvar_splice_stat" and len(path) == 1)
                  or (cache == "alphamissense_tranx_domain_map" and not path)
                  or (cache == "clinvar_gene_stat" and
                      (not path or (len(path) == 2 and path[0] in ("domains", "transcripts"))))
                  or (cache == "pm1_regions_pkl" and not path and
                      top_key in ("hcseeker_hotspots", "rmc_constrained")))
        if record:
            pairs.sort(key=lambda pair: serialized(pair[0]))
        return [kind, pairs]
    if kind in ("list", "tuple", "set", "frozenset"):
        items = [normalize(cache, top_key, item, path + (index,)) for index, item in enumerate(value[1])]
        if kind in ("set", "frozenset"):
            items.sort(key=serialized)
        return [kind, items]
    return value


def reference_value(cache: str, key: object, value: object) -> list:
    # Schema10 explicitly represents a float NaN as a None code, not as a numeric score.
    if cache == "clingen_map" and type(value) is float and math.isnan(value):
        value = None
    return normalize(cache, key, encode(value))


def differences(expected: object, actual: object, path: tuple = ()):
    """Yield every unequal primitive or container span without truncating the ledger."""
    if type(expected) is not type(actual):
        yield path, expected, actual
    elif isinstance(expected, list):
        if len(expected) != len(actual):
            yield path + ("length",), len(expected), len(actual)
        for index in range(max(len(expected), len(actual))):
            if index >= len(expected):
                yield path + (index,), {"absent": True}, actual[index]
            elif index >= len(actual):
                yield path + (index,), expected[index], {"absent": True}
            else:
                yield from differences(expected[index], actual[index], path + (index,))
    elif expected != actual:
        yield path, expected, actual


def sha256(path: Path) -> str:
    with path.open("rb") as stream:
        return hashlib.file_digest(stream, "sha256").hexdigest()


def compare(args: argparse.Namespace) -> dict:
    started = time.monotonic()
    identities = {name: {"path": str(path), "bytes": path.stat().st_size, "sha256": sha256(path)}
                  for name, path in (("pickle", args.pickle), ("archive", args.archive), ("binary", args.binary),
                                     ("comparator", Path(__file__).resolve()))}
    with args.pickle.open("rb") as stream:
        compressed = stream.read(2) == b"\x1f\x8b"
    opener = gzip.open if compressed else open
    with opener(args.pickle, "rb") as stream:
        original = pickle.load(stream)
        if stream.read(1):
            raise ValueError("reference has trailing bytes after the pickle")
    if args.cache == "all_intolerant_domains":
        if type(original) is not set:
            raise TypeError("all_intolerant_domains must be a set")
        expected_by_key = {None: original}
    else:
        if type(original) is not dict:
            raise TypeError("reference root must be a dict")
        expected_by_key = original
    seen = set()
    last_key = None
    compared = 0
    mismatch_count = 0
    expected_digest, actual_digest = hashlib.sha256(), hashlib.sha256()
    args.report.parent.mkdir(parents=True, exist_ok=True)
    args.mismatches.parent.mkdir(parents=True, exist_ok=True)
    errors_path = args.report.with_suffix(".stderr.txt")
    with args.mismatches.open("x") as mismatches, errors_path.open("x") as errors:
        def record(key, reason, expected=None, actual=None, path=()):
            nonlocal mismatch_count
            mismatch_count += 1
            json.dump({"key": key, "reason": reason, "path": path,
                       "expected": expected, "actual": actual}, mismatches, ensure_ascii=False)
            mismatches.write("\n")
        command = [str(args.binary), "cache-export", args.cache, "--input", str(args.archive), "--output", "-"]
        process = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=errors, text=True, encoding="utf-8")
        try:
            assert process.stdout is not None
            for line in process.stdout:
                key, value = json.loads(line)
                decoded = key_value(key)
                if decoded in seen:
                    record(key, "duplicate top-level key")
                if compared and decoded is not None and last_key is not None and decoded <= last_key:
                    record(key, "top-level keys not strictly sorted")
                last_key = decoded
                seen.add(decoded)
                actual = normalize(args.cache, decoded, value)
                actual_digest.update((serialized([key, actual]) + "\n").encode())
                if decoded not in expected_by_key:
                    record(key, "extra top-level key")
                else:
                    expected = reference_value(args.cache, decoded, expected_by_key[decoded])
                    expected_digest.update((serialized([key, expected]) + "\n").encode())
                    for path, left, right in differences(expected, actual):
                        record(key, "value mismatch", left, right, path)
                compared += 1
        except BaseException:
            process.kill()
            process.wait()
            raise
        finally:
            if process.stdout is not None:
                process.stdout.close()
        returncode = process.wait()
        for key in expected_by_key:
            if key not in seen:
                record(encode(key), "missing top-level key")
        if returncode:
            record(None, "native export failed", 0, returncode)
    stable = all(sha256(Path(item["path"])) == item["sha256"] for item in identities.values())
    result = {"cache": args.cache, "identities": identities, "inputs_unchanged": stable,
              "python": sys.version,
              "reference_entries": len(expected_by_key), "exported_entries": compared,
              "unique_exported_keys": len(seen), "mismatches": mismatch_count,
              "reference_values_sha256": expected_digest.hexdigest(), "native_values_sha256": actual_digest.hexdigest(),
              "native_returncode": returncode, "wall_seconds": time.monotonic() - started,
              "reference_peak_rss_kib": resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
              "native_peak_rss_kib": resource.getrusage(resource.RUSAGE_CHILDREN).ru_maxrss,
              "passed": stable and not mismatch_count and not returncode and compared == len(expected_by_key)}
    with args.report.open("x") as report:
        json.dump(result, report, indent=2)
        report.write("\n")
    print(json.dumps({k: v for k, v in result.items() if k != "identities"}), flush=True)
    return result


class ComparatorTests(unittest.TestCase):
    def test_float_bits_and_scalar_types(self):
        nan = struct.unpack(">d", bytes.fromhex("7ff8000000000042"))[0]
        self.assertEqual(encode(nan), ["float", "7ff8000000000042"])
        self.assertNotEqual(encode(-0.0), encode(0.0))
        self.assertNotEqual(encode(True), encode(1))
        self.assertNotEqual(encode(1), encode(1.0))

    def test_containers_and_missingness(self):
        self.assertNotEqual(encode([]), encode(()))
        self.assertNotEqual(encode(set()), encode(frozenset()))
        self.assertNotEqual(encode({"x": None}), encode({}))
        self.assertNotEqual(encode(["x", "x"]), encode(["x"]))

    def test_only_fixed_record_order_is_ignored(self):
        first = {"pp_domains": [], "exon_domains": {"2": [], "1": []}}
        fields_reordered = {"exon_domains": {"2": [], "1": []}, "pp_domains": []}
        exons_reordered = {"pp_domains": [], "exon_domains": {"1": [], "2": []}}
        normalized = lambda x: reference_value("alphamissense_tranx_domain_map", "t", x)
        self.assertEqual(normalized(first), normalized(fields_reordered))
        self.assertNotEqual(normalized(first), normalized(exons_reordered))

    def test_spliceai_name_order_and_sequence_order_are_protected(self):
        a = [{"pos": 1, "splice_ai": {"a": "1", "b": "2"}}]
        b = [{"splice_ai": {"a": "1", "b": "2"}, "pos": 1}]
        c = [{"pos": 1, "splice_ai": {"b": "2", "a": "1"}}]
        normalized = lambda x: reference_value("clinvar_splice_stat", "t", x)
        self.assertEqual(normalized(a), normalized(b))
        self.assertNotEqual(normalized(a), normalized(c))
        self.assertNotEqual(normalized(a + c), normalized(c + a))

    def test_pm1_lookup_keys_only_are_sorted(self):
        a, b = {"b": {2}, "a": {1}}, {"a": {1}, "b": {2}}
        self.assertEqual(reference_value("pm1_regions_pkl", "hcseeker_hotspots", a), reference_value("pm1_regions_pkl", "hcseeker_hotspots", b))
        self.assertNotEqual(reference_value("pm1_regions_pkl", "metadata", a), reference_value("pm1_regions_pkl", "metadata", b))

    def test_every_difference_is_enumerated(self):
        found = list(differences([1, 2, 3], [4, 5]))
        self.assertEqual(len(found), 4)  # length, two changed leaves, one missing item
        self.assertEqual(reference_value("clingen_map", "v", float("nan")), ["none"])
        self.assertNotEqual(reference_value("clinvar_patho_af_stat", "g", float("nan")), ["none"])


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    commands = parser.add_subparsers(dest="command", required=True)
    p = commands.add_parser("compare")
    for name in ("binary", "pickle", "archive", "report", "mismatches"):
        p.add_argument("--" + name, type=Path, required=True)
    p.add_argument("--cache", required=True)
    commands.add_parser("selftest")
    args = parser.parse_args()
    if args.command == "selftest":
        suite = unittest.defaultTestLoader.loadTestsFromTestCase(ComparatorTests)
        raise SystemExit(not unittest.TextTestRunner(verbosity=2).run(suite).wasSuccessful())
    raise SystemExit(0 if compare(args)["passed"] else 1)


if __name__ == "__main__":
    main()
