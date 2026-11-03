#!/usr/bin/env python3
"""Compare the whole native DAS step with the original analyze_domain_data.

check --binary PRIVA --source-dir SOURCE --work-root NEW
oracle runs analyze_domain_data in its own interpreter under PYTHONHASHSEED
equal to the native --python-hash-seed. Only the two rescue candidate iterators
are replaced, by the approved domain-sorted-rescue-v1 reference
(das_sorted_reference.py); every other procedure is the unchanged source.

Each case builds one synthetic AlphaMissense VEP VCF that both the Python
statistics script (final pickle) and the native builder (typed Parquet) collect,
an InterPro mapping pickle with its native archive and, when the case needs
them, native domain metadata, a coldspot table and a ClinVar VCF. Approved
differences are normalised: tolerance rows are keyed by domain (Python writes
worker completion order, native query order) with every row compared as text;
the pickled sets are compared as sorted lists. The filter table is compared as
text. Transcript mappings, plots and metadata caches are not compared.
"""
from __future__ import annotations
import argparse
import ast
import difflib
import gzip
import hashlib
import importlib
import importlib.util
import json
import logging
import os
from pathlib import Path
import pickle
import random
import shutil
import subprocess
import sys
from types import SimpleNamespace

sys.dont_write_bytecode = True
SORTED_PIN = "1a2b4efc312c41c6fb20afa3c1d90566da7cd0772d2f261c2f7f09d094fae834"
ASSEMBLY = "hg19"
THREADS = 2
FDR = "0.01"
OVERLAP = "0.6"
HIGH, INTOLERANT, LOW = (0.86, 0.99), (0.93, 0.99), (0.05, 0.3)
GO_MF_MIDDLE, NON_GO_MF_MIDDLE = (0.6, 0.75), (0.76, 0.84)
SPOTS = ("type\tgene\taa_start_pos\taa_end_pos\n"
         "coldspot\tSYM_COLDRM\t1\t30\n"
         "coldspot\tSYM_COLDKEEP\t10\t100\n"
         "hotspot\tSYM_INTDOMAIN\t1\t100\n"
         "coldspot\tSYM_ABSENT\t1\t5\n")
NORMALIZATION = {
    "database_name_mappings": {"vep_to_interpro": {"Pfam": "PFAM", "SMART": "SMART"},
                               "interpro_to_vep": {"PFAM": "Pfam", "SMART": "SMART"}},
    "accession_id_normalization": {"vep_to_interpro_key": {"__default__": {"rule_type": "identity"}},
                                   "interpro_key_to_vep_id": {"__default__": {"rule_type": "identity"}}},
}
CURATED = ("IPR_ID\tname\ttype\tshort_name\tsource\tgo_terms\tnotes\tMolecular_Function_GO_Terms\tFunctionality_Assessment\n"
           "IPR000001\tname\tDomain\tshort\tsynthetic\t.\t.\tGO:0003674\tFunctional\n")
COUNT_FIELDS = ["domains", "references", "queries", "ks_intolerant", "coldspot_removed", "benign_removed",
                "common_removed", "family_removed", "rescued_go_mf", "rescued_non_go_mf", "final_intolerant"]
FILTERS = ["coldspot", "benign", "common", "family"]


def load(path, name):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def domain(label, kind=None, mf=False, bounds=HIGH, count=7, residues=None, length=200, afs=None, scores=None):
    residues = list(residues) if residues is not None else list(range(10, 10 + count))
    if scores is None:
        rng = random.Random("das-assembly:" + label)
        scores = [f"{rng.uniform(*bounds):.4f}" for _ in residues]
    if len(scores) != len(residues):
        raise ValueError("one score per residue is required")
    return dict(label=label, kind=kind, mf=mf, residues=residues, length=length, scores=scores,
                afs=afs or ["0.0"] * len(residues), symbol="SYM_" + label)


def full_domains():
    references = [domain(f"REF{i}", "Domain", True) for i in range(6)] + [domain("REFFAMILY", "Family", True)]
    bimodal = [f"{0.10 + 0.01 * i:.4f}" for i in range(9)] + [f"{0.95 + 0.004 * i:.4f}" for i in range(10)]
    queries = [
        domain("TOLORPHAN", bounds=LOW),
        domain("TOLDOMAIN", "Domain", bounds=LOW, count=10),
        domain("INTDOMAIN", "Domain", bounds=INTOLERANT),
        domain("INTORPHAN", bounds=INTOLERANT),
        domain("COLDRM", "Homologous_superfamily", bounds=INTOLERANT, residues=range(1, 8)),
        domain("COLDKEEP", "Family", bounds=INTOLERANT, residues=range(1, 13)),
        domain("BENRM", "Domain", bounds=INTOLERANT),
        domain("BENKEEP", "Domain", bounds=INTOLERANT),
        domain("COMMONRM", "Domain", bounds=INTOLERANT, afs=["0.2"] + ["0.0"] * 6),
        domain("FAMRM", "Family", True, INTOLERANT, residues=[1, 60, 120, 180, 240, 300, 360, 400], length=400),
        domain("FAMKEEP", "Family", bounds=INTOLERANT, length=400),
        domain("GOMFREPEAT", "Repeat", True, GO_MF_MIDDLE, residues=[1, 50, 100, 150, 200, 250, 300], length=300),
        domain("NONGOMF", "Domain", bounds=NON_GO_MF_MIDDLE),
        domain("NONGOMFORPHAN", bounds=NON_GO_MF_MIDDLE),
        domain("BIMODAL", "Repeat", residues=range(1, 20), scores=bimodal),
    ]
    return references + queries


def type_domains():
    references = [domain(f"T{i:03d}", "Domain", True, count=6) for i in range(105)]
    references += [domain(f"TFAMILY{i}", "Family", True, count=6, length=100) for i in range(3)]
    queries = [
        domain("TINTDOMAIN", "Domain", bounds=INTOLERANT),
        domain("TTOLDOMAIN", "Domain", bounds=LOW),
        domain("TNONGOMF", "Domain", bounds=NON_GO_MF_MIDDLE),
        domain("TINTFAMILY", "Family", bounds=INTOLERANT),
        domain("TNONGOMFORPHAN", bounds=NON_GO_MF_MIDDLE),
        domain("TINTORPHAN", bounds=INTOLERANT),
    ]
    return references + queries


def cases():
    full = full_domains()
    return {
        "full": dict(domains=full, metadata=True, spots=True, clinvar=True, seed=1,
                     expect=dict(tolerant=1, coldspot_removed=1, benign_removed=1, common_removed=1,
                                 family_removed=1, rescued_go_mf=1, rescued_non_go_mf=1)),
        "no_metadata": dict(domains=full, seed=0, expect=dict(tolerant=1, common_removed=1)),
        "type_composite": dict(domains=type_domains(), metadata=True, seed=42,
                               expect=dict(composites=2, tolerant=1, rescued_non_go_mf=1)),
        "no_queries": dict(domains=[domain(f"ONLYREF{i}", "Domain", True) for i in range(6)], metadata=True,
                           seed=4294967295, fail=True),
    }


def build_vcf(factory, folder, domains):
    rows, position = [], 100
    for item in domains:
        label = item["label"]
        for residue, score, af in zip(item["residues"], item["scores"], item["afs"]):
            csq = factory.annotation(Gene="ENSG_" + label, Feature="ENST_" + label, EXON="1/1",
                                     Protein_position=f"{residue}/{item['length']}", DOMAINS="Pfam:PF" + label,
                                     SYMBOL=item["symbol"])
            rows.append(factory.row(pos=position, score=score, af=af, variant=f"A{residue}T", csq=csq))
            position += 1
    return factory.make_case(folder, {"rows": rows})


def build_mapping(normalizer, domains):
    mapping = {}
    for item in domains:
        if item["kind"] is None:
            continue
        key = normalizer.normalize_vep_to_interpro("Pfam", "PF" + item["label"])
        if key is None:
            raise RuntimeError("fixture normalization profile changed")
        term = (("GO:0003674", "molecular_function", "molecular function") if item["mf"]
                else ("GO:0008150", "biological_process", "biological process"))
        mapping[key] = [("IPR_" + item["label"], item["kind"], "short", "name", [term], "abstract")]
    return mapping


def build_clinvar(benign, helper, folder):
    rows, position = [], 100

    def add(label, **change):
        nonlocal position
        csq = benign.annotation(Gene="ENSG_" + label, DOMAINS="Pfam:PF" + label, Feature="ENST_" + label)
        rows.append(benign.row(pos=position, csq=csq, **change))
        position += 1
    for _ in range(3):
        add("BENRM")
    for _ in range(2):
        add("BENKEEP")
    add("BENKEEP", sig="Likely_benign")
    add("BENKEEP", review="criteria_provided,_single_submitter")
    return helper.make_case(folder, {"fields": benign.FIELDS, "rows": rows})


def source_tree(args):
    data = args.source_dir / "data/InterPro"
    normalization, curated = data / "InterPro_domain_norm_map.json", data / "curated_InterPro_func_domains.tsv.gz"
    if normalization.is_file() and curated.is_file():
        return args.source_dir / "scripts", normalization, "source"
    tree = args.work_root / "source_tree"
    shutil.copytree(args.source_dir / "scripts", tree / "scripts", ignore=shutil.ignore_patterns("__pycache__"))
    data = tree / "data/InterPro"
    data.mkdir(parents=True)
    (data / "InterPro_domain_norm_map.json").write_text(json.dumps(NORMALIZATION, indent=2) + "\n")
    (data / "curated_InterPro_func_domains.tsv.gz").write_bytes(gzip.compress(CURATED.encode(), mtime=0))
    return tree / "scripts", data / "InterPro_domain_norm_map.json", "synthetic"


def run(command, folder, name, env=None, timeout=1800):
    with (folder / (name + ".stdout.log")).open("xb") as out, (folder / (name + ".stderr.log")).open("xb") as err:
        return subprocess.run(list(map(str, command)), stdout=out, stderr=err, env=env, timeout=timeout).returncode


def oracle(args):
    sys.path.insert(0, str(args.scripts))
    das = importlib.import_module("am_pick_intolerant_domains")
    reference = load(args.sorted_reference, "sorted_rescue_reference")
    source = Path(das.__file__)
    rendered = reference.render_source(source.read_bytes())
    if hashlib.sha256(rendered).hexdigest() != SORTED_PIN:
        raise RuntimeError("sorted rescue reference differs from the approved profile")
    functions = [node for node in ast.parse(rendered).body
                 if isinstance(node, ast.FunctionDef) and node.name in reference.FUNCTIONS]
    if len(functions) != len(reference.FUNCTIONS):
        raise RuntimeError("sorted rescue procedures absent")
    exec(compile(ast.fix_missing_locations(ast.Module(body=functions, type_ignores=[])), str(source), "exec"), vars(das))
    counts = dict(composites=0)

    def observe(name, update):
        original = getattr(das, name)

        def observed(*positional, **named):
            value = original(*positional, **named)
            update(value)
            return value
        setattr(das, name, observed)
    observe("visualize_domain_distribution", lambda value: counts.update(domains=len(value)))
    observe("build_fixed_equal_weight_composite", lambda value: counts.update(composites=counts["composites"] + 1))
    observe("rescue_bimodal_go_mf_domains", lambda value: counts.update(rescued_go_mf=len(value[0])))
    observe("rescue_non_gomf_domains", lambda value: counts.update(rescued_non_go_mf=len(value[0])))

    def optional(value):
        return None if value is None else str(value)
    completed = das.analyze_domain_data(str(args.statistics), str(args.output_dir), threads=args.threads,
                                        fdr_threshold=args.fdr_threshold, assembly=ASSEMBLY,
                                        am_vcf=optional(args.am_vcf), hcseeker_spots_tsv=optional(args.spots),
                                        clinvar_vcf=optional(args.clinvar_vcf),
                                        coldspot_overlap_threshold=args.coldspot_overlap_threshold,
                                        interpro_mapping_pickle=str(args.mapping))
    result = dict(completed=completed is True, counts=counts)
    if completed is True:
        def read(name):
            with (args.output_dir / f"{name}.{ASSEMBLY}.pkl").open("rb") as stream:
                return pickle.load(stream)
        result["all_intolerant_domains"] = sorted(read("all_intolerant_domains"))
        result["curated_functional_references"] = sorted(read("curated_functional_references"))
        result["excluded_domains"] = [[name, sorted(values)] for name, values in read("excluded_domains").items()]
    args.result.write_text(json.dumps(result) + "\n")
    if completed is not True:
        raise SystemExit(3)


def table(path):
    if not path.is_file():
        return None
    text = path.read_bytes().decode()
    lines = text.split("\n")
    terminated = lines[-1] == ""
    lines = lines[:-1] if terminated else lines
    rows, duplicates = {}, []
    for line in lines[1:]:
        key = line.split("\t", 1)[0]
        if key in rows:
            duplicates.append(key)
        rows[key] = line
    return dict(header=lines[0] if lines else None, rows=rows, duplicates=duplicates, terminated=terminated)


def exported_domains(path):
    lines = [json.loads(line) for line in path.read_text().splitlines() if line]
    if len(lines) != 1 or lines[0][0] != ["none"] or lines[0][1][0] != "set":
        raise ValueError("unexpected all_intolerant_domains export layout")
    values = []
    for item in lines[0][1][1]:
        if item[0] != "str":
            raise ValueError("non-text intolerant domain")
        values.append(item[1])
    if len(set(values)) != len(values):
        raise ValueError("duplicate intolerant domain")
    return sorted(values)


def compare(args, folder, python, problems, matched):
    python_root, native_root = folder / "python_output", folder / "native_output"

    def tsv(root, name):
        return root / f"{name}.{ASSEMBLY}.tsv"
    left, right = (table(tsv(root, "domain_tolerance_analysis")) for root in (python_root, native_root))
    if left is None or right is None:
        matched["domain_tolerance_analysis"] = False
        problems.append("tolerance table missing")
    else:
        same = (left["header"] == right["header"] and left["terminated"] and right["terminated"]
                and not left["duplicates"] and not right["duplicates"] and left["rows"] == right["rows"])
        matched["domain_tolerance_analysis"] = same
        if not same:
            problems.append("tolerance tables differ")
            differing = {key: dict(python=value, native=right["rows"].get(key))
                         for key, value in left["rows"].items() if right["rows"].get(key) != value}
            (folder / "tolerance_difference.json").write_text(json.dumps(dict(
                header=dict(python=left["header"], native=right["header"]),
                terminated=dict(python=left["terminated"], native=right["terminated"]),
                duplicates=dict(python=left["duplicates"], native=right["duplicates"]),
                only_native=sorted(set(right["rows"]) - set(left["rows"])), differing=differing), indent=1) + "\n")
    filters = [tsv(root, "domain_filter_results") for root in (python_root, native_root)]
    if not all(path.is_file() for path in filters):
        matched["domain_filter_results"] = False
        problems.append("filter table missing")
    else:
        texts = [path.read_bytes().decode() for path in filters]
        matched["domain_filter_results"] = texts[0] == texts[1]
        if texts[0] != texts[1]:
            problems.append("filter tables differ")
            (folder / "filter_results.diff").write_text("\n".join(difflib.unified_diff(
                texts[0].split("\n"), texts[1].split("\n"), "python", "native", lineterm="")) + "\n")
    differences = {}
    archive = native_root / f"all_intolerant_domains.{ASSEMBLY}.rkyv"
    native_sets = {}
    if not archive.is_file() or run([args.binary, "cache-export", "all_intolerant_domains", "--input", archive,
                                     "--output", "-"], folder, "native_export"):
        problems.append("native intolerant archive missing or unreadable")
    else:
        try:
            native_sets["all_intolerant_domains"] = exported_domains(folder / "native_export.stdout.log")
        except (ValueError, IndexError, TypeError) as error:
            problems.append("native intolerant export: " + str(error))
    for name in ["curated_functional_references", "excluded_domains"]:
        path = native_root / f"{name}.{ASSEMBLY}.json"
        if not path.is_file():
            problems.append(name + " missing")
            continue
        value = json.loads(path.read_text())
        if name == "excluded_domains":
            if list(value) != [key for key in FILTERS if key in value]:
                problems.append("excluded filter order differs from coldspot, benign, common, family")
            value = [[key, members] for key, members in value.items()]
            arrays = [members for _, members in value]
        else:
            arrays = [value]
        if any(members != sorted(set(members)) for members in arrays):
            problems.append(name + " arrays are not sorted unique text")
        native_sets[name] = value
    for name in ["all_intolerant_domains", "curated_functional_references", "excluded_domains"]:
        if name not in native_sets:
            matched[name] = False
            continue
        matched[name] = native_sets[name] == python[name]
        if not matched[name]:
            problems.append(name + " differs")
            differences[name] = dict(python=python[name], native=native_sets[name])
    if differences:
        (folder / "set_difference.json").write_text(json.dumps(differences, indent=1) + "\n")
    return left, filters


def expected_counts(python, tolerance, filter_path):
    excluded = dict(python["excluded_domains"])
    filter_rows = [line for line in filter_path.read_text().split("\n")[1:] if line] if filter_path.is_file() else []
    counts = dict(domains=python["counts"].get("domains"), references=len(python["curated_functional_references"]),
                  queries=len(tolerance["rows"]) if tolerance else None, ks_intolerant=len(filter_rows))
    for name in FILTERS:
        counts[name + "_removed"] = len(excluded[name]) if name in excluded else None
    counts.update(rescued_go_mf=python["counts"].get("rescued_go_mf"),
                  rescued_non_go_mf=python["counts"].get("rescued_non_go_mf"),
                  final_intolerant=len(python["all_intolerant_domains"]))
    return counts


def native_counts(folder):
    lines = [line for line in (folder / "native_analysis.stdout.log").read_text().splitlines() if line.strip()]
    return json.loads(lines[-1]) if lines else None


def run_case(args, tools, name, spec):
    folder = args.work_root / name
    folder.mkdir()
    record, problems, matched = dict(case=name, hash_seed=spec["seed"]), [], {}
    am = folder / "alphamissense"
    am.mkdir()
    vcf = build_vcf(tools.factory, am, spec["domains"])
    statistics_pickle, statistics = folder / "statistics.pkl", folder / "statistics.parquet"
    mapping_pickle, mapping = folder / "mapping.pkl", folder / "mapping.rkyv"
    steps = dict(
        python_statistics=run([sys.executable, tools.tests / "alpha_domain_compare.py", "oracle", "--procedure", "statistics",
                               "--source", tools.scripts / "stat_protein_domain_amscores.py", "--vcf", vcf,
                               "--output", statistics_pickle, "--threads", "1"], folder, "python_statistics"),
        native_statistics=run([args.binary, "alphamissense-domain", "build", "--vcf", vcf, "--output", statistics,
                               "--threads", "1"], folder, "native_statistics"))
    with mapping_pickle.open("xb") as stream:
        pickle.dump(build_mapping(tools.normalizer, spec["domains"]), stream, protocol=4)
    steps["native_mapping"] = run([args.binary, "cache-convert", "interpro_mapping_pickle", "--pickle", mapping_pickle,
                                   "--output", mapping], folder, "native_mapping")
    python_options, native_options = [], []
    if spec.get("metadata"):
        metadata = folder / "metadata.parquet"
        steps["native_metadata"] = run([args.binary, "domain-metadata", "build", "--vcf", vcf, "--output", metadata],
                                       folder, "native_metadata")
        python_options += ["--am-vcf", vcf]
        native_options += ["--metadata", metadata]
    if spec.get("spots"):
        spots = folder / "spots.tsv"
        spots.write_text(SPOTS)
        python_options += ["--spots", spots]
        native_options += ["--spots", spots]
    if spec.get("clinvar"):
        clinvar_folder = folder / "clinvar"
        clinvar_folder.mkdir()
        clinvar = build_clinvar(tools.benign, tools.clinvar, clinvar_folder)
        python_options += ["--clinvar-vcf", clinvar]
        native_options += ["--clinvar-vcf", clinvar]
    record["setup_exits"] = steps
    if any(steps.values()):
        record.update(passed=False, problems=["input preparation failed"])
        return record
    import pyarrow.parquet as pq
    builder = pq.ParquetFile(statistics).metadata.metadata[b"priva.builder_source_sha256"].decode()
    seed = str(spec["seed"])
    environment = dict(os.environ, PYTHONHASHSEED=seed, MPLBACKEND="Agg", PYTHONDONTWRITEBYTECODE="1")
    python_result = folder / "python_result.json"
    python_exit = run([sys.executable, Path(__file__).resolve(), "oracle", "--scripts", tools.scripts,
                       "--sorted-reference", tools.tests / "das_sorted_reference.py", "--statistics", statistics_pickle,
                       "--mapping", mapping_pickle, "--output-dir", folder / "python_output", "--result", python_result,
                       "--threads", THREADS, "--fdr-threshold", FDR, "--coldspot-overlap-threshold", OVERLAP,
                       *python_options], folder, "python_analysis", environment)
    native_exit = run([args.binary, "das-domains", "analyze", "--statistics", statistics, "--source-builder-sha256", builder,
                       "--interpro", mapping, "--normalization-json", tools.normalization, "--python-hash-seed", seed,
                       "--threads", THREADS, "--assembly", ASSEMBLY, "--output-dir", folder / "native_output",
                       "--fdr-threshold", FDR, "--coldspot-overlap-threshold", OVERLAP, *native_options],
                      folder, "native_analysis")
    record.update(python_exit=python_exit, native_exit=native_exit)
    if spec.get("fail"):
        if python_exit != 3:
            problems.append("Python analysis did not end with its no-reference/no-query return")
        if native_exit == 0:
            problems.append("native analysis succeeded where Python returned None")
        published = sorted(path.name for path in (folder / "native_output").glob("*")) if (folder / "native_output").is_dir() else []
        if published:
            problems.append("native analysis published outputs after failing")
        if (folder / "python_output" / f"all_intolerant_domains.{ASSEMBLY}.pkl").exists():
            problems.append("Python wrote an intolerant set despite returning None")
        record.update(native_published=published, passed=not problems, problems=problems)
        return record
    if python_exit or native_exit:
        problems.append("analysis failed: python exit %d, native exit %d" % (python_exit, native_exit))
        record.update(passed=False, problems=problems)
        return record
    python = json.loads(python_result.read_text())
    tolerance, filters = compare(args, folder, python, problems, matched)
    expected = expected_counts(python, tolerance, filters[0])
    actual = native_counts(folder)
    record.update(python_counts=expected, native_counts=actual, filters_run=[key for key, _ in python["excluded_domains"]])
    matched["counts"] = actual is not None and list(actual) == COUNT_FIELDS and actual == expected
    if not matched["counts"]:
        problems.append("count line differs")
        (folder / "count_difference.json").write_text(json.dumps(dict(python=expected, native=actual), indent=1) + "\n")
    observed = dict(expected, composites=python["counts"].get("composites"),
                    tolerant=(expected["queries"] or 0) - expected["ks_intolerant"])
    record["unexercised"] = [key for key, minimum in spec.get("expect", {}).items() if (observed.get(key) or 0) < minimum]
    record["composites_built"] = observed["composites"]
    record.update(matched=matched, passed=not problems, problems=problems)
    return record


def check(args):
    logging.disable(logging.CRITICAL)
    args.work_root.mkdir(parents=True, exist_ok=False)
    tests = args.source_dir / "rust/tests"
    scripts, normalization, origin = source_tree(args)
    lookup = load(tests / "das_domain_lookup_compare.py", "assembly_lookup")
    tools = SimpleNamespace(
        tests=tests, scripts=scripts, normalization=normalization,
        factory=load(tests / "alpha_domain_compare.py", "assembly_alphamissense_fixture"),
        clinvar=load(tests / "clinvar_aa_splice_compare.py", "assembly_clinvar_fixture"),
        benign=load(tests / "das_benign_compare.py", "assembly_benign_fixture"),
        normalizer=lookup.normalizer(scripts / "protein_domain_mapping.py", json.loads(normalization.read_text())))
    results = []
    for name, spec in cases().items():
        record = run_case(args, tools, name, spec)
        results.append(record)
        print(json.dumps(record), flush=True)
    summary = dict(passed=all(record["passed"] for record in results), cases=len(results), interpro_data=origin,
                   rescue_profile="domain-sorted-rescue-v1", results=results)
    (args.work_root / "complete.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps({key: value for key, value in summary.items() if key != "results"}), flush=True)
    if not summary["passed"]:
        raise SystemExit(1)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    commands = parser.add_subparsers(dest="command", required=True)
    test = commands.add_parser("check")
    for name in ("binary", "source-dir", "work-root"):
        test.add_argument("--" + name, type=Path, required=True)
    original = commands.add_parser("oracle")
    for name in ("scripts", "sorted-reference", "statistics", "mapping", "output-dir", "result"):
        original.add_argument("--" + name, type=Path, required=True)
    for name in ("am-vcf", "spots", "clinvar-vcf"):
        original.add_argument("--" + name, type=Path)
    original.add_argument("--threads", type=int, required=True)
    original.add_argument("--fdr-threshold", type=float, required=True)
    original.add_argument("--coldspot-overlap-threshold", type=float, required=True)
    args = parser.parse_args()
    globals()[args.command](args)


if __name__ == "__main__":
    main()
