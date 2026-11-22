#!/usr/bin/env python3
"""Compare `priva spliceai2-prepare` with its unchanged Python source.

check --binary PRIVA --repo CHECKOUT --work-root NEW_DIRECTORY [--manifest REAL_MANIFEST
      --vep-annotation GFF3 [--subset-chrom chr21]] [--only CASE ...]
oracle --producer SCRIPT --exception JSON -- ARGUMENTS (internal reference launcher)

Each case builds two identical fixture folders. The Python arm runs
CHECKOUT/scripts/spliceai2_to_vep_vcf.py through the oracle launcher, which
records the class and text of the exception its main() reports through
parser.exit; the native arm runs `PRIVA spliceai2-prepare`.

- gene-context cases use generated GTF/GFF3 and scored-region fixtures that
  reach every check and the VEP-release renaming rules.
- prepare failure cases use generated manifests and stop before conversion.
- With --manifest, `real_gene_context` maps the real scored regions against
  the real GENCODE48 GTF and VEP annotation, and `real_subset_prepare` runs a
  --validation-subset preparation of the SNV and indel shards of one
  chromosome with the real inputs (read only). The Python arm builds its
  standalone converter with rustc as the script does; the native arm uses
  its own spliceai2-vcf.

Compared exactly: exit status; stdout lines that report verification,
projection or a JSON result; the exception category and message (case folder
written as <ARM>); every written file (context TSVs and VCF/index bytes by
SHA-256); and receipts, except `converter` and `driver_sha256`, which name
the program that converted by design. Argument errors exit 2 in both.
Writes WORK_ROOT/complete.json and exits 1 on any mismatch.
"""
from __future__ import annotations

import argparse
import gzip
import hashlib
import importlib.util
import json
import os
from pathlib import Path
import re
import shutil
import subprocess
import sys

sys.dont_write_bytecode = True
PRODUCER_SHA256 = "548f29e85a4be78ab96c5afe95009a56e6967d47df8b78842ff60156b37ad813"
DESCRIPTION = "##description: evidence-based annotation of the human genome (GRCh38), version 48 (Ensembl 114)\n"
VEP_DESCRIPTION = "##description: evidence-based annotation of the human genome (GRCh38), version 50 (Ensembl 116)\n"
REGION_HEADER = "gene_id\tchrom\tgene_start\tgene_end\tstrand\n"
REVISION = "cd1b8b82bcef9c97ac3b3b9a24f0de3a5bbf77f9"
KEPT = ("Verifying ", "Projecting ", "{")


def load(path, name):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def write_json(path, value):
    Path(path).write_text(json.dumps(value, sort_keys=True, indent=2) + "\n")


def oracle(args):
    """Run the unchanged main() and record the exception it reports."""
    module = load(args.producer, "original_spliceai2_to_vep_vcf")
    record = {}
    original_exit = module.argparse.ArgumentParser.exit

    def observed_exit(parser, status=0, message=None):
        kind, value, _ = sys.exc_info()
        if status == 1 and kind is not None:
            record.update(category=kind.__name__, message=str(value))
            write_json(args.exception, record)
        return original_exit(parser, status, message)

    module.argparse.ArgumentParser.exit = observed_exit
    sys.argv = [str(args.producer), *args.arguments]
    try:
        module.main()
    except SystemExit:
        raise
    except BaseException as exc:
        record.update(category=type(exc).__name__, message=str(exc), uncaught=True)
        write_json(args.exception, record)
        raise
    raise SystemExit(0)


# ---------------------------------------------------------------- fixtures

def gene(gene_id, name, start, end, strand="+", chrom="chr1", kind="protein_coding", extra=""):
    attributes = f'gene_id "{gene_id}"; gene_type "{kind}"; gene_name "{name}"; level 2;{extra}'
    return f"{chrom}\tHAVANA\tgene\t{start}\t{end}\t.\t{strand}\t.\t{attributes}\n"


def exon(gene_id, start, end, chrom="chr1"):
    return f'{chrom}\tHAVANA\texon\t{start}\t{end}\t.\t+\t.\tgene_id "{gene_id}"; transcript_id "T";\n'


BASE_GENES = [gene("ENSG00000000001.5", "ALPHA", 100, 900), exon("ENSG00000000001.5", 100, 200),
              gene("ENSG00000000002.3", "BETA", 1000, 2000, "-"),
              gene("ENSG00000000003.1", "GAMMA", 50, 80, "+", "chr2"),
              gene("ENSG00000000004.2", "DELTA", 5000, 6000, "+", "chrX"),
              gene("ENSG00000000004.2_PAR_Y", "DELTA", 5000, 6000, "+", "chrY"),
              gene("ENSG00000000009.1", "NONCODE", 10, 20, kind="lncRNA")]


def gtf(genes=None, description=DESCRIPTION):
    return description + "".join(BASE_GENES if genes is None else genes)


def regions(rows):
    return REGION_HEADER + "".join("\t".join(map(str, row)) + "\n" for row in rows)


GOOD_REGIONS = [("ENSG00000000001.5", "chr1", 100, 900, "+"), ("ENSG00000000002.3", "chr1", 1000, 2000, "-"),
                ("ENSG00000000003.1", "chr2", 50, 80, "+")]


def vep_gtf(rows, description=VEP_DESCRIPTION):
    return description + "".join(gene(*row) for row in rows)


RENAMES = [("ENSG00000000001.6", "ALPHA", 100, 900), ("ENSG00000000002.4", "BETA2", 1000, 2000, "-"),
           ("ENSG00000000003.2", "GAMMA", 50, 80, "-", "chr2")]


def context_case(region_rows=None, genes=None, description=DESCRIPTION, vep=None, release="116",
                 vep_name="vep.gtf", gz=False, extra_args=()):
    files = {"regions.tsv": regions(GOOD_REGIONS if region_rows is None else region_rows),
             "gencode.gtf.gz" if gz else "gencode.gtf": ("gzip", gtf(genes, description)) if gz else gtf(genes, description)}
    args = ["gene-context", "--regions", "{arm}/regions.tsv",
            "--gtf", "{arm}/gencode.gtf" + (".gz" if gz else ""), "--output", "{arm}/context.tsv"]
    if vep is not None:
        files[vep_name] = vep
        args += ["--vep-annotation", "{arm}/" + vep_name, "--ensembl-release", release]
    return {"files": files, "args": [*args, *extra_args]}


def catalogue(real):
    cases = {
        "argument_unknown_command": {"args": ["convert"]},
        "argument_threads_text": {"args": ["prepare", "--manifest", "/m", "--assembly", "hg38", "--output-dir", "/o",
                                           "--receipt", "/o/r.json", "--threads", "two"]},
        "context_exact": context_case(),
        "context_gzip_gtf": context_case(gz=True),
        "context_zero_based": context_case([("ENSG00000000001.5", "chr1", 99, 900, "+"),
                                            ("ENSG00000000002.3", "chr1", 999, 2000, "-")]),
        "context_mixed_coordinates": context_case([("ENSG00000000001.5", "chr1", 99, 900, "+"),
                                                   ("ENSG00000000002.3", "chr1", 1000, 2000, "-")]),
        "context_stable_ids": context_case([("ENSG00000000001", "chr1", 100, 900, "+"),
                                            ("ENSG00000000004_PAR_Y", "chrY", 5000, 6000, "+"),
                                            ("ENSG00000000004.2", "chrX", 5000, 6000, "+")]),
        "context_ambiguous_stable_id": context_case([("ENSG00000000001", "chr1", 100, 900, "+")],
                                                    BASE_GENES + [gene("ENSG00000000001.6", "ALPHA2", 100, 900)]),
        "context_two_sources_one_gene": context_case([("ENSG00000000001", "chr1", 100, 900, "+"),
                                                      ("ENSG00000000001.5", "chr1", 100, 900, "+")]),
        "context_duplicate_region": context_case([GOOD_REGIONS[0], GOOD_REGIONS[0]]),
        "context_unknown_gene": context_case([("ENSG00000000077.1", "chr1", 1, 2, "+")]),
        "context_strand_mismatch": context_case([("ENSG00000000001.5", "chr1", 100, 900, "-")]),
        "context_not_protein_coding": context_case([("ENSG00000000009.1", "chr1", 10, 20, "+")]),
        "context_bounds_differ": context_case([("ENSG00000000001.5", "chr1", 100, 901, "+")]),
        "context_unsafe_symbol": context_case(None, BASE_GENES[:1] + [gene("ENSG00000000002.3", "BE TA", 1000, 2000, "-"),
                                                                    BASE_GENES[3]]),
        "context_opposite_strand_symbol": context_case(
            [GOOD_REGIONS[0], ("ENSG00000000008.1", "chr1", 3000, 4000, "-")],
            BASE_GENES + [gene("ENSG00000000008.1", "ALPHA", 3000, 4000, "-")]),
        "context_same_strand_symbol": context_case(
            [GOOD_REGIONS[0], ("ENSG00000000008.1", "chr1", 3000, 4000, "+")],
            BASE_GENES + [gene("ENSG00000000008.1", "ALPHA", 3000, 4000, "+")]),
        "context_empty_regions": context_case([]),
        "context_wrong_schema": {**context_case(), "files": {**context_case()["files"], "regions.tsv": "gene\tchrom\n"}},
        "context_malformed_region": {**context_case(), "files": {**context_case()["files"],
                                                                 "regions.tsv": REGION_HEADER + "ENSG00000000001.5\tchr1\t100\n"}},
        "context_other_gencode": context_case(description=DESCRIPTION.replace("version 48", "version 47")),
        "context_duplicate_gtf_id": context_case(genes=BASE_GENES + [BASE_GENES[0]]),
        "context_missing_gene_type": context_case(genes=[gene("ENSG00000000001.5", "ALPHA", 100, 900).replace(
            'gene_type "protein_coding"; ', "")]),
        "context_duplicate_gene_name": context_case(genes=[gene("ENSG00000000001.5", "ALPHA", 100, 900,
                                                                extra=' gene_name "ALPHA";')]),
        "context_short_gtf_row": context_case(genes=["chr1\tHAVANA\tgene\t1\t2\n"]),
        "context_bad_coordinate": context_case(genes=[gene("ENSG00000000001.5", "ALPHA", "1x", 900)]),
        "context_release_names_gtf": context_case(vep=vep_gtf(RENAMES)),
        "context_release_names_gff3": context_case(
            vep=VEP_DESCRIPTION.replace("##", "#") + "".join(
                f"{c}\tENSEMBL\tgene\t{s}\t{e}\t.\t{st}\t.\tID=gene:{i};gene_id={i};gene_name={n}\n"
                for i, n, s, e, st, c in [("ENSG00000000001.6", "AL%50HA", 100, 900, "+", "chr1"),
                                          ("ENSG00000000002.4", "BETA%2C2", 1000, 2000, "-", "chr1")]),
            vep_name="vep.gff3"),
        "context_release_name_conflicts": context_case(
            [*GOOD_REGIONS, ("ENSG00000000008.1", "chr1", 3000, 4000, "+")],
            BASE_GENES + [gene("ENSG00000000008.1", "THETA", 3000, 4000)],
            vep=vep_gtf([("ENSG00000000001.6", "SHARED", 100, 900), ("ENSG00000000008.2", "SHARED", 3000, 4000),
                         ("ENSG00000000002.4", "ALPHA", 1000, 2000, "-")])),
        "context_release_not_found": context_case(vep=vep_gtf([("ENSG00000000005.1", "OTHER", 1, 2)])),
        "context_release_mismatch": context_case(vep=vep_gtf(RENAMES), release="115"),
        "context_release_wrong_format": context_case(vep=vep_gtf(RENAMES), vep_name="vep.txt"),
        "context_release_without_number": {**context_case(), "args": [*context_case()["args"],
                                                                      "--vep-annotation", "{arm}/regions.tsv"]},
        "context_output_exists": {**context_case(), "files": {**context_case()["files"], "context.tsv": "kept\n"}},
    }
    small = {"files": {"m.json": json.dumps({"schema_version": 1})}}
    prepare = ["prepare", "--manifest", "{arm}/m.json", "--assembly", "GRCh38", "--output-dir", "{arm}/out",
               "--receipt", "{arm}/out/adapter.ready.json"]

    def manifest(**changes):
        value = {"schema_version": 1, "status": "downloads_complete", "confirmed_by": "fixture",
                 "resource": "SpliceAI2", "version": "2.0", "assembly": "GRCh38", "source_revision": REVISION,
                 "gencode_gtf": {"release": "GENCODE48"}, "vep_gene_annotation": {},
                 "files": [{"kind": "snv", "chrom": "chr21", "path": "{arm}/SpliceAI2_GRCh38_snv_chr21.tsv.gz"},
                           {"kind": "indel", "chrom": "chr21", "path": "{arm}/SpliceAI2_GRCh38_indel_chr21.tsv.gz"}]}
        value.update(changes)
        return {"files": {"m.json": json.dumps(value), "SpliceAI2_GRCh38_snv_chr21.tsv.gz": "x",
                          "SpliceAI2_GRCh38_indel_chr21.tsv.gz": "x"}, "args": prepare}

    cases.update({
        "prepare_hg19": {**small, "args": [*prepare[:4], "hg19", *prepare[5:]]},
        "prepare_zero_threads": {**small, "args": [*prepare, "--threads", "0"]},
        "prepare_relative_manifest": {**small, "args": ["prepare", "--manifest", "m.json", *prepare[3:]]},
        "prepare_partial_manifest": {"files": {"m.partial.json": "{}"}, "args": [prepare[0], prepare[1], "{arm}/m.partial.json", *prepare[3:]]},
        "prepare_downloading_manifest": manifest(status="downloading"),
        "prepare_wrong_revision": manifest(source_revision="other"),
        "prepare_other_gencode": manifest(gencode_gtf={"release": "GENCODE47"}),
        "prepare_without_vep_annotation": manifest(vep_gene_annotation=None),
        "prepare_full_set_required": manifest(),
        "prepare_subset_one_kind": {**manifest(files=[{"kind": "snv", "chrom": "chr21", "path": "{arm}/SpliceAI2_GRCh38_snv_chr21.tsv.gz"}]),
                                    "args": [*prepare, "--validation-subset"]},
        "prepare_unpublished_shard": manifest(files=[{"kind": "indel", "chrom": "chrY", "path": "/x"}]),
        "prepare_duplicate_shard": manifest(files=[{"kind": "snv", "chrom": "chr21", "path": "/x"}] * 2),
        "prepare_shard_without_kind": manifest(files=[{"chrom": "chr21", "path": "/x"}]),
        "prepare_wrong_source_name": {**manifest(files=[{"kind": "snv", "chrom": "chr21", "path": "{arm}/chr21.tsv.gz"},
                                                        {"kind": "indel", "chrom": "chr21", "path": "{arm}/SpliceAI2_GRCh38_indel_chr21.tsv.gz"}]),
                                      "args": [*prepare, "--validation-subset"]},
        "prepare_regions_size_mismatch": {**manifest(scored_regions={"path": "{arm}/regions.tsv", "bytes": 5, "sha256": "0" * 64}),
                                          "args": [*prepare, "--validation-subset"]},
    })
    cases["prepare_regions_size_mismatch"]["files"]["regions.tsv"] = regions(GOOD_REGIONS)
    if real:
        cases["real_gene_context"] = {"real": "context"}
        cases["real_subset_prepare"] = {"real": "prepare"}
    return cases


def real_case(kind, arm, real_manifest, vep_annotation, chrom):
    data = json.loads(real_manifest.read_text())
    annotation = {"path": str(vep_annotation), "bytes": vep_annotation.stat().st_size,
                  "sha256": hashlib.sha256(vep_annotation.read_bytes()).hexdigest(), "release": "GENCODE50"}
    if kind == "context":
        return ["gene-context", "--regions", data["scored_regions"]["path"], "--gtf", data["gencode_gtf"]["path"],
                "--output", str(arm / "context.tsv"), "--vep-annotation", str(vep_annotation),
                "--ensembl-release", "116"]
    data["vep_gene_annotation"] = annotation
    data["files"] = [item for item in data["files"] if item["chrom"] == chrom]
    (arm / "subset.json").write_text(json.dumps(data, indent=2) + "\n")
    return ["prepare", "--manifest", str(arm / "subset.json"), "--assembly", "GRCh38",
            "--output-dir", str(arm / "out"), "--receipt", str(arm / "out/adapter.ready.json"),
            "--threads", "4", "--validation-subset"]


# ---------------------------------------------------------------- runs

def materialize(arm, case):
    arm.mkdir(parents=True)
    for relative, spec in case.get("files", {}).items():
        path = arm / relative
        path.parent.mkdir(parents=True, exist_ok=True)
        if isinstance(spec, tuple) and spec[0] == "gzip":
            path.write_bytes(gzip.compress(spec[1].replace("{arm}", str(arm)).encode(), mtime=0))
        else:
            path.write_text(spec.replace("{arm}", str(arm)))


def normalize(text, arm):
    if not isinstance(text, str):
        return text
    for spelling in sorted({str(arm), os.path.realpath(arm)}, key=len, reverse=True):
        text = text.replace(spelling, "<ARM>")
    return text


def receipt_view(data, arm):
    body = json.loads(data)
    body.pop("converter", None)
    body.pop("driver_sha256", None)
    return json.loads(normalize(json.dumps(body), arm))


def tree(root):
    state = {}
    for folder, directories, names in os.walk(root, followlinks=False):
        directories[:] = [d for d in directories if d != ".tools"]
        for name in sorted(directories + names):
            path = Path(folder) / name
            relative = normalize(str(path.relative_to(root)), root)
            if path.is_symlink():
                state[relative] = ["link", normalize(os.readlink(path), root)]
            elif path.is_dir():
                state[relative] = ["directory", oct(path.stat().st_mode & 0o7777)]
            elif name.endswith(".ready.json"):
                state[relative] = ["receipt", receipt_view(path.read_bytes(), root), oct(path.stat().st_mode & 0o7777)]
            elif path.stat().st_size <= 65536 and not name.endswith((".gz", ".tbi")):
                state[relative] = ["text", normalize(path.read_text(errors="replace"), root), oct(path.stat().st_mode & 0o7777)]
            else:
                digest = hashlib.sha256()
                with path.open("rb") as handle:
                    for block in iter(lambda: handle.read(1 << 23), b""):
                        digest.update(block)
                state[relative] = ["bytes", digest.hexdigest(), oct(path.stat().st_mode & 0o7777)]
    return state


def run(command, folder, label, environment, cwd):
    with (folder / (label + ".stdout")).open("wb") as out, (folder / (label + ".stderr")).open("wb") as err:
        process = subprocess.run(list(map(str, command)), cwd=cwd, stdout=out, stderr=err,
                                 stdin=subprocess.DEVNULL, env=environment, timeout=6 * 3600)
    return process.returncode


def outcome(folder, label, status, arm):
    stdout = (folder / (label + ".stdout")).read_bytes().decode("utf-8", "replace")
    stderr = (folder / (label + ".stderr")).read_bytes().decode("utf-8", "replace")
    kept = [normalize(line, arm) for line in stdout.splitlines() if line.startswith(KEPT)]
    result = {"exit": status, "stdout": kept, "tree": tree(arm)}
    if label == "python":
        recorded = folder / "exception.json"
        value = json.loads(recorded.read_text()) if recorded.exists() else {}
        if value:
            result.update(category=value.get("category"), message=value.get("message"))
    else:
        matched = re.search(r"(?:\A|\n)priva: ([A-Za-z_][A-Za-z_0-9]*): ((?:(?!\npriva: ).)*)\n\Z", stderr, re.S)
        if matched:
            result.update(category=matched[1], message=matched[2])
        elif status not in (0, 2):
            result.update(category=None, message=stderr[-4000:])
    if "message" in result:
        result["message"] = normalize(result["message"], arm)
    return result


def compare(python, native, before):
    if python["exit"] == 2 or native["exit"] == 2:
        return python["exit"] == native["exit"] == 2 and python["tree"] == native["tree"] == before
    left, right = dict(python), dict(native)
    left.pop("uncaught", None)
    return left == right


def check(args):
    producer = args.repo / "scripts/spliceai2_to_vep_vcf.py"
    if hashlib.sha256(producer.read_bytes()).hexdigest() != PRODUCER_SHA256:
        raise SystemExit("spliceai2_to_vep_vcf.py differs from the qualified source; update the pin")
    cases = catalogue(args.manifest is not None)
    unknown = set(args.only or ()) - cases.keys()
    if unknown:
        raise SystemExit("Unknown cases: " + ", ".join(sorted(unknown)))
    args.work_root.mkdir(parents=True, exist_ok=False)
    rows = {}
    for name, case in cases.items():
        if args.only and name not in args.only:
            continue
        folder = args.work_root / name
        folder.mkdir()
        results, befores = {}, {}
        for label in ("python", "native"):
            arm = (folder / label).absolute()
            materialize(arm, case)
            if "real" in case:
                arguments = real_case(case["real"], arm, args.manifest, args.vep_annotation, args.subset_chrom)
            else:
                arguments = [value.replace("{arm}", str(arm)) for value in case["args"]]
            befores[label] = tree(arm)
            environment = {key: value for key, value in os.environ.items() if not key.startswith("PRIVA_")}
            environment["PYTHONDONTWRITEBYTECODE"] = "1"
            if label == "python":
                command = [sys.executable, Path(__file__).resolve(), "oracle", "--producer", producer,
                           "--exception", folder / "exception.json", "--", *arguments]
            else:
                command = [args.binary, "spliceai2-prepare", *arguments]
            status = run(command, folder, label, environment, arm)
            results[label] = outcome(folder, label, status, arm)
        passed = befores["python"] == befores["native"] and \
            compare(results["python"], results["native"], befores["native"])
        rows[name] = {"passed": passed, "python": results["python"], "native": results["native"]}
        write_json(folder / "comparison.json", rows[name])
        if "real" in case and passed:
            for label in ("python", "native"):
                shutil.rmtree(folder / label)
    failed = sorted(name for name, row in rows.items() if not row["passed"])
    report = {
        "producer": str(producer), "producer_sha256": PRODUCER_SHA256,
        "binary": str(args.binary),
        "binary_sha256": hashlib.sha256(Path(args.binary).read_bytes()).hexdigest(),
        "cases": len(rows), "python_successes": sum(row["python"]["exit"] == 0 for row in rows.values()),
        "failed": failed, "rows": rows,
    }
    write_json(args.work_root / "complete.json", report)
    print(json.dumps({key: report[key] for key in ("cases", "python_successes", "failed")}))
    if failed:
        raise SystemExit(1)


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    commands = parser.add_subparsers(dest="procedure", required=True)
    p = commands.add_parser("check")
    p.add_argument("--binary", type=Path, required=True)
    p.add_argument("--repo", type=Path, required=True)
    p.add_argument("--work-root", type=Path, required=True)
    p.add_argument("--manifest", type=Path)
    p.add_argument("--vep-annotation", type=Path)
    p.add_argument("--subset-chrom", default="chr21")
    p.add_argument("--only", action="append")
    p.set_defaults(function=check)
    p = commands.add_parser("oracle")
    p.add_argument("--producer", type=Path, required=True)
    p.add_argument("--exception", type=Path, required=True)
    p.add_argument("arguments", nargs=argparse.REMAINDER)
    p.set_defaults(function=oracle)
    args = parser.parse_args()
    if args.procedure == "oracle" and args.arguments[:1] == ["--"]:
        args.arguments = args.arguments[1:]
    args.function(args)


if __name__ == "__main__":
    main()
