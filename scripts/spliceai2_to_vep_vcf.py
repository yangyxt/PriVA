#!/usr/bin/env python3
"""Prepare a verified, annotation-only SpliceAI2 view; never activate ACMG policy.

Rust performs the row conversion. This deployment helper verifies immutable
inputs, qualifies gene context against the same GENCODE release, and manages
the mature gzip/bgzip/tabix tools. Original 84-column tables remain authoritative.
"""
import argparse
import csv
import fcntl
import gzip
import hashlib
import json
import os
from pathlib import Path
import re
import shutil
import subprocess
import tempfile

REVISION = "cd1b8b82bcef9c97ac3b3b9a24f0de3a5bbf77f9"
PLUGIN_SHA256 = "500a534c403bcbcc15753718768f0a1f352cc6b93f9119e5eba25745f7205b9a"
CHROMS = [f"chr{i}" for i in range(1, 23)] + ["chrX", "chrY"]
PUBLISHED = {(kind, chrom) for kind in ("snv", "indel") for chrom in CHROMS
             if (kind, chrom) != ("indel", "chrY")}
REGION_COLUMNS = ["gene_id", "chrom", "gene_start", "gene_end", "strand"]
CONTEXT_COLUMNS = REGION_COLUMNS + ["gene_name"]


def require(condition, message):
    if not condition:
        raise ValueError(message)


def digest(path):
    h = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(8 * 1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def identity(path):
    s = path.stat()
    return (s.st_dev, s.st_ino, s.st_size, s.st_mtime_ns, s.st_ctime_ns)


def complete_file(path):
    path = Path(path)
    require(path.is_absolute(), f"An absolute input path is required: {path}")
    require(not any(".partial" in part.lower() for part in path.parts), f"Incomplete input: {path}")
    resolved = path.resolve(strict=True)
    require(not any(".partial" in part.lower() for part in resolved.parts), f"Incomplete input target: {path}")
    require(resolved.is_file() and resolved.stat().st_size > 0, f"Missing/empty input: {path}")
    return resolved


def verified_file(item, identities):
    path = complete_file(item["path"])
    require(type(item.get("bytes")) is int and item["bytes"] > 0, "Expected positive input bytes")
    require(re.fullmatch(r"[0-9a-f]{64}", item.get("sha256", "")), "Expected pinned SHA256")
    before = identity(path)
    require(before[2] == item["bytes"], f"Input size mismatch: {path}")
    print(f"Verifying {path.name}", flush=True)
    require(digest(path) == item["sha256"], f"Input SHA256 mismatch: {path}")
    require(identity(path) == before, f"Input changed during verification: {path}")
    identities[path] = before
    return path


def text_stream(path):
    return gzip.open(path, "rt", encoding="utf-8", newline="") if path.suffix == ".gz" else path.open(encoding="utf-8", newline="")


def gene_context(regions, gtf, destination):
    """Require exact gene IDs and bounds, not symbol aliases or guessed intervals."""
    genes, description = {}, []
    with text_stream(gtf) as stream:
        for number, line in enumerate(stream, 1):
            if line.startswith("#"):
                description.append(line)
                continue
            fields = line.rstrip("\r\n").split("\t")
            require(len(fields) == 9, f"Invalid GTF row {number}")
            if fields[2] != "gene":
                continue
            attributes = re.findall(r'(\w+) "([^"\n]*)";', fields[8])
            for key in ("gene_id", "gene_name", "gene_type"):
                require(sum(k == key for k, _ in attributes) == 1, f"Missing/duplicate GTF {key} at {number}")
            attributes = dict(attributes)
            gene_id, symbol = attributes.get("gene_id"), attributes.get("gene_name")
            require(gene_id and symbol, f"GTF gene lacks ID/name at {number}")
            require(gene_id not in genes, f"Duplicate GTF gene ID: {gene_id}")
            genes[gene_id] = (fields[0], int(fields[3]), int(fields[4]), fields[6], symbol,
                              attributes.get("gene_type"))
    description = "".join(description)
    require("GRCh38" in description and "version 48" in description,
            "Gene mapping requires the original GRCh38 GENCODE version 48 GTF")
    rows, conventions, seen, symbols = [], set(), set(), {}
    with text_stream(regions) as stream:
        reader = csv.DictReader(stream, delimiter="\t")
        require(reader.fieldnames == REGION_COLUMNS, "Unexpected publisher scored-region schema")
        for row in reader:
            require(None not in row and all(value is not None for value in row.values()), "Malformed scored-region row")
            gene_id = row["gene_id"]
            require(gene_id not in seen, f"Duplicate scored-region gene: {gene_id}")
            seen.add(gene_id)
            require(gene_id in genes, f"No exact GENCODE48 gene ID: {gene_id}; version stripping/aliases are not allowed")
            chrom, start, end, strand, name, gene_type = genes[gene_id]
            require(chrom == row["chrom"] and strand == row["strand"], f"Gene contig/strand mismatch: {gene_id}")
            require(gene_type == "protein_coding", f"Scored region is not a GENCODE48 protein-coding gene: {gene_id}")
            raw_start, raw_end = int(row["gene_start"]), int(row["gene_end"])
            require(raw_end == end and raw_start in (start, start - 1), f"Scored-region/GTF bounds differ: {gene_id}")
            conventions.add("1-based-closed" if raw_start == start else "0-based-half-open")
            require(chrom in CHROMS and strand in ("+", "-") and 0 < start <= end, f"Unsupported gene context: {gene_id}")
            require(not re.search(r"[\s|;,=]", name), f"Gene symbol cannot be represented safely in the plugin: {gene_id}")
            require((chrom, name) not in symbols or symbols[chrom, name] == gene_id,
                    f"Multiple gene IDs share a symbol on {chrom}: {name}; plugin matching would be ambiguous")
            symbols[chrom, name] = gene_id
            rows.append((gene_id, chrom, start, end, strand, name))
    require(rows and len(conventions) == 1, "Empty or mixed-coordinate scored-region table")
    rows.sort(key=lambda row: (CHROMS.index(row[1]), row[2], row[3], row[0]))
    with destination.open("x", newline="") as stream:
        writer = csv.writer(stream, delimiter="\t", lineterminator="\n")
        writer.writerow(CONTEXT_COLUMNS)
        writer.writerows(rows)
    return {"genes": len(rows), "source_coordinates": conventions.pop(),
            "context_coordinates": "1-based-closed", "gene_id_matching": "exact",
            "join": "source variant position and model strand within verified scored-gene intervals",
            "scope": "Adapter-derived gene association, not a publisher gene-conditioned prediction. Predicted sites may lie outside the associated gene; no effects are filtered."}


def terminate(processes):
    for process in processes:
        if process.poll() is None:
            process.terminate()
    for process in processes:
        try:
            process.wait(timeout=5)
        except subprocess.TimeoutExpired:
            process.kill()
            process.wait()
        for stream in (process.stdin, process.stdout, process.stderr):
            if stream is not None and not stream.closed:
                try:
                    stream.close()
                except BrokenPipeError:
                    pass  # Preserve the original child/process failure.


def convert_kind(kind, entries, context, target, converter, threads, tools):
    """Sequential chromosome writers share one compressor; no Python row copying."""
    processes, statistics = [], []
    try:
        with target.open("xb") as output:
            compressor = subprocess.Popen([tools["bgzip"], "-c", "-@", str(threads)], stdin=subprocess.PIPE, stdout=output)
            processes.append(compressor)
            for index, item in enumerate(entries):
                print(f"Projecting {kind} {item['chrom']}", flush=True)
                source = subprocess.Popen([tools["gzip"], "-dc", str(item["resolved_path"])], stdout=subprocess.PIPE)
                processes.append(source)
                command = [str(converter), "spliceai2-vcf", "--genes", str(context), "--chrom", item["chrom"],
                           "--kind", kind, "--input", "-", "--output", "-"]
                if index:
                    command.append("--no-header")
                projection = subprocess.Popen(command, stdin=source.stdout, stdout=compressor.stdin, stderr=subprocess.PIPE)
                processes.append(projection)
                source.stdout.close()
                _, diagnostic = projection.communicate()
                require(projection.returncode == 0, "Rust projection failed: " + diagnostic.decode(errors="replace")[-4096:])
                require(source.wait() == 0, f"Source decompression failed: {item['path']}")
                report = json.loads(diagnostic)
                require(report["module"] == "spliceai2-vcf" and report["chrom"] == item["chrom"] and report["kind"] == kind,
                        "Invalid Rust conversion report")
                require(report["input_rows"] > 0 and report["output_records"] >= report["input_rows"], "Empty/incomplete projection")
                statistics.append(report)
            compressor.stdin.close()
            require(compressor.wait() == 0, "BGZF compression failed")
        subprocess.run([tools["tabix"], "-p", "vcf", str(target)], check=True)
        contigs = subprocess.check_output([tools["tabix"], "-l", str(target)], text=True).splitlines()
        require(contigs == [item["chrom"] for item in entries], "Indexed contigs differ from converted inputs")
        # Probe the first observed record, not an arbitrary possibly empty interval.
        with gzip.open(target, "rt") as stream:
            first = next(line for line in stream if not line.startswith("#")).rstrip("\n")
        fields = first.split("\t")
        query = subprocess.check_output([tools["tabix"], str(target), f"{fields[0]}:{fields[1]}-{fields[1]}"], text=True)
        require(first in query.splitlines(), "Nonempty VCF/index query did not recover the first record")
        return statistics
    finally:
        terminate(processes)


def prepare(args):
    require(args.assembly in ("GRCh38", "hg38"), "SpliceAI2 route supports GRCh38/hg38 only; hg19 remains v1")
    require(args.threads > 0, "threads must be positive")
    manifest_path = complete_file(args.manifest)
    manifest_before = identity(manifest_path)
    manifest = json.loads(manifest_path.read_text())
    manifest_sha = digest(manifest_path)
    require(type(manifest.get("schema_version")) is int and manifest["schema_version"] == 1 and
            manifest.get("status") == "downloads_complete" and isinstance(manifest.get("confirmed_by"), str) and manifest["confirmed_by"].strip(),
            "A confirmed, complete source manifest is required; downloading inventories are not inputs")
    require(manifest.get("resource") == "SpliceAI2" and manifest.get("version") == "2.0" and
            manifest.get("assembly") == "GRCh38" and manifest.get("source_revision") == REVISION,
            "Wrong source version/assembly/revision")
    require(manifest.get("gencode_gtf", {}).get("release") == "GENCODE48", "Pinned GENCODE48 GTF is required")
    files = manifest.get("files", [])
    pairs = [(item["kind"], item["chrom"]) for item in files]
    require(len(set(pairs)) == len(pairs) and set(pairs) <= PUBLISHED and pairs, "Duplicate/unsupported source shards")
    require(set(pairs) == PUBLISHED or args.validation_subset, "Full cache requires all 47 published shards; chrY indels are not published")
    require(not args.validation_subset or {kind for kind, _ in pairs} == {"snv", "indel"}, "Subset validation requires both variant types")
    driver = Path(__file__).resolve()
    identities = {manifest_path: manifest_before, driver: identity(driver)}
    driver_sha = digest(driver)
    for item in files:
        expected = f"SpliceAI2_GRCh38_{item['kind']}_{item['chrom']}.tsv.gz"
        require(Path(item["path"]).name == expected, f"Unexpected source name: {item['path']}")
        complete_file(item["path"])  # Cheap completeness checks before large reads.
    regions = verified_file(manifest["scored_regions"], identities)
    gtf = verified_file(manifest["gencode_gtf"], identities)
    plugin = verified_file(manifest["vep_spliceai_plugin"], identities)
    require(manifest["vep_spliceai_plugin"]["sha256"] == PLUGIN_SHA256, "Unqualified official plugin revision")
    converter_name = shutil.which(str(args.priva_bin))
    require(converter_name, "Build the Rust priva executable and supply --priva-bin")
    converter = complete_file(Path(converter_name).absolute())
    identities[converter] = identity(converter)
    converter_sha = digest(converter)
    tool_identities, tools = {}, {}
    for name in ("gzip", "bgzip", "tabix"):
        executable = shutil.which(name)
        require(executable, f"Missing mature command tool: {name}")
        path = Path(executable).resolve(strict=True)
        tools[name] = str(path)
        identities[path] = identity(path)
        version = subprocess.check_output([str(path), "--version"], stderr=subprocess.STDOUT, text=True).splitlines()[0]
        tool_identities[name] = {"path": str(path), "sha256": digest(path), "version": version}
    folder = args.output_dir.absolute()
    receipt_path = args.receipt.absolute()
    for path in (folder, receipt_path):
        require(not any(".partial" in part.lower() or part == ".." for part in path.parts), "Incomplete/noncanonical output path")
        require(not any(parent.is_symlink() for parent in (path, *path.parents)), "Output paths must not traverse symlinks")
    require(receipt_path.parent.resolve() == folder.resolve(), "Receipt must be directly inside the output directory")
    folder.mkdir(parents=True, exist_ok=True)
    final_paths = [folder / f"SpliceAI2_GRCh38.{kind}.vep.vcf.gz{suffix}" for kind in ("snv", "indel") for suffix in ("", ".tbi")]
    context_name = "SpliceAI2_GRCh38.gene_context.tsv"
    final_paths += [folder / context_name, receipt_path]
    with (folder / ".install.lock").open("a") as lock:
        fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        require(not any(path.exists() or path.is_symlink() for path in final_paths), "Refusing to replace existing prepared outputs")
        with tempfile.TemporaryDirectory(prefix=".building-", dir=folder) as temporary:
            stage = Path(temporary)
            context = stage / context_name
            mapping = gene_context(regions, gtf, context)
            for item in files:
                item["resolved_path"] = verified_file(item, identities)
            outputs, statistics = {}, {}
            for kind in ("snv", "indel"):
                entries = sorted((item for item in files if item["kind"] == kind), key=lambda item: CHROMS.index(item["chrom"]))
                target = stage / f"SpliceAI2_GRCh38.{kind}.vep.vcf.gz"
                statistics[kind] = convert_kind(kind, entries, context, target, converter, args.threads, tools)
                index = Path(str(target) + ".tbi")
                outputs[kind] = {"path": str(folder / target.name), "index_path": str(folder / index.name),
                                 "bytes": target.stat().st_size, "sha256": digest(target), "index_sha256": digest(index)}
            for path, before in identities.items():
                require(identity(path) == before, f"Input changed while preparing cache: {path}")
            receipt = {
                "schema_version": 1, "resource": "SpliceAI2", "version": "2.0", "assembly": "GRCh38",
                "status": "validated_subset_annotation_only" if args.validation_subset else "prepared_annotation_only",
                "coverage": "validation_subset" if args.validation_subset else "complete_published",
                "scientific_qualified": False, "source_revision": REVISION,
                "source_manifest_sha256": manifest_sha, "outputs": outputs, "statistics": statistics,
                "source_files": [{k: v for k, v in item.items() if k != "resolved_path"} for item in files],
                "scored_regions": manifest["scored_regions"], "gencode_gtf": manifest["gencode_gtf"],
                "gene_context": {"path": str(folder / context.name), "sha256": digest(context), **mapping},
                "vep_spliceai_plugin": {"path": str(plugin), "sha256": PLUGIN_SHA256},
                "converter": {"path": str(converter), "sha256": converter_sha}, "tools": tool_identities,
                "driver_sha256": driver_sha, "slurm_job_id": os.environ.get("SLURM_JOB_ID"),
                "thresholds": {"high_recall": 0.1, "balanced": 0.25, "high_precision": 0.5},
                "projection": "Unchanged top0 site delta/distance strings only, one allele/gene per VCF record. No cutoff filtering, score recomputation, masking, normalization or strand flipping. Original 84-column tables retain all top3/junction/reference/alternate/summary values.",
                "coverage_limits": "GENCODE48 protein-coding scored regions; population-observed indels; no published chrY indel shard. Missing coverage is not zero score.",
                "license": {"source": f"https://huggingface.co/datasets/illumina-ai/SpliceAI2-data/tree/{REVISION}",
                            "note": "Original and derived score resources retain Illumina SpliceAI2 terms; no redistribution or public service is enabled by this adapter."},
                "activation": "Preparation only. Legacy PriVA ACMG/ClinVar/SpliceVault and CADD v1 feature contracts are not qualified for substitution."}
            staged_receipt = stage / receipt_path.name
            require(staged_receipt not in (context, *[stage / p.name for p in final_paths[:-2]]), "Receipt name collides with an output")
            staged_receipt.write_text(json.dumps(receipt, indent=2) + "\n")
            # Receipt is published last; incomplete sets never have a readiness receipt.
            for path in final_paths[:-1]:
                (stage / path.name).replace(path)
            staged_receipt.replace(receipt_path)
    print(json.dumps({"status": receipt["status"], "receipt": str(receipt_path)}), flush=True)
    return receipt


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest", type=Path, required=True)
    parser.add_argument("--assembly", required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--receipt", type=Path, required=True)
    parser.add_argument("--priva-bin", default="priva")
    parser.add_argument("--threads", type=int, default=1)
    parser.add_argument("--validation-subset", action="store_true", help="Explicit bounded test; never produces an installable full-cache receipt")
    args = parser.parse_args()
    try:
        prepare(args)
    except (ValueError, OSError, subprocess.SubprocessError, KeyError) as error:
        parser.exit(1, f"SpliceAI2 preparation failed: {error}\n")


if __name__ == "__main__":
    main()
