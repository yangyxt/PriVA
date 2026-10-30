#!/usr/bin/env python3
"""Prepare a verified, annotation-only SpliceAI2 view; never activate ACMG policy.

The standalone installation Rust tool performs the row conversion. This helper verifies immutable
inputs, qualifies gene context against the same GENCODE release, and manages
the mature gzip/bgzip/tabix tools. Original 84-column tables remain authoritative.
The official plugin attaches a score only when its gene name equals VEP's, so each
gene is also labelled with its name in the GENCODE release of VEP's Ensembl release.
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
import sys
import tempfile
from urllib.parse import unquote

REVISION = "cd1b8b82bcef9c97ac3b3b9a24f0de3a5bbf77f9"
QUALIFIED_PLUGINS = {
    "500a534c403bcbcc15753718768f0a1f352cc6b93f9119e5eba25745f7205b9a": {
        "release": "113", "allele_matching": "Ensembl get_matched_variant_alleles normalization"},
    "d08ce7d8c2e6b3229cef638beff145e465e224f0e3fcd43fd49d501ebc2740b1": {
        "release": "116", "allele_matching": "Exact position/REF/ALT after VF indel conversion and VF-strand orientation"},
}
CHROMS = [f"chr{i}" for i in range(1, 23)] + ["chrX", "chrY"]
PUBLISHED = {(kind, chrom) for kind in ("snv", "indel") for chrom in CHROMS
             if (kind, chrom) != ("indel", "chrY")}
REGION_COLUMNS = ["gene_id", "chrom", "gene_start", "gene_end", "strand"]
CONTEXT_COLUMNS = REGION_COLUMNS + ["gene_name"]
MIN_RUST_VERSION = (1, 70, 0)  # is_some_and/is_ok_and; no edition-2024 requirement.
ISOLATED_RUST_VERSION = "1.90.0"


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


def compiler_identity(executable):
    """Return the actual compatible compiler identity, or None."""
    if not executable:
        return None
    executable = shutil.which(str(executable))
    if not executable:
        return None
    # Keep the rustc filename: rustup and other multicall shims dispatch by argv[0].
    path = Path(executable).absolute()
    try:
        version = subprocess.check_output([str(path), "--version"], text=True, stderr=subprocess.STDOUT).strip()
        sysroot = subprocess.check_output([str(path), "--print", "sysroot"], text=True, stderr=subprocess.STDOUT).strip()
    except (OSError, subprocess.SubprocessError):
        return None
    match = re.match(r"rustc ([0-9]+)\.([0-9]+)\.([0-9]+)(?:\s|$)", version)
    if not match or tuple(map(int, match.groups())) < MIN_RUST_VERSION:
        return None
    compiler = Path(sysroot) / "bin/rustc"
    if not compiler.is_file():
        return None
    compiler = compiler.resolve(strict=True)
    return {"path": str(path), "version": version, "executable": str(compiler),
            "sha256": digest(compiler), "invocation_sha256": digest(path), "sysroot": sysroot}


def installation_compiler(tools_dir, rustc=None):
    """Use a compatible host compiler or acquire a private pinned Conda compiler."""
    compiler = compiler_identity(rustc or "rustc")
    if compiler:
        return compiler
    require(rustc is None, "The explicit --rustc must be an executable Rust 1.70+ compiler")
    prefix = tools_dir / ("rust-" + ISOLATED_RUST_VERSION)
    with (tools_dir / ".compiler.lock").open("a") as lock:
        fcntl.flock(lock, fcntl.LOCK_EX)
        compiler = compiler_identity(prefix / "bin/rustc")
        if compiler:
            require(compiler["version"].split()[1] == ISOLATED_RUST_VERSION,
                    "The private compiler prefix does not match pinned Rust " + ISOLATED_RUST_VERSION)
            return compiler
        require(not prefix.exists(), "The isolated Rust prefix is incomplete; explicit recovery is required: " + str(prefix))
        candidates = [os.environ.get("MAMBA_EXE"), os.environ.get("CONDA_EXE"), "mamba", "micromamba", "conda"]
        manager = next((resolved for candidate in candidates if candidate
                        for resolved in [shutil.which(candidate)] if resolved), None)
        require(manager, "SpliceAI2 installation needs a compatible rustc or the existing Conda/Mamba manager")
        print("Installing isolated SpliceAI2 build compiler: " + str(prefix), flush=True)
        package_env = dict(os.environ, CONDA_PKGS_DIRS=str(tools_dir / ".conda-pkgs"))
        subprocess.run([manager, "create", "--yes", "--prefix", str(prefix), "--override-channels",
                        "--channel", "conda-forge", "rust=" + ISOLATED_RUST_VERSION],
                       env=package_env, check=True)
        compiler = compiler_identity(prefix / "bin/rustc")
        require(compiler and compiler["version"].split()[1] == ISOLATED_RUST_VERSION,
                "The isolated compiler did not provide pinned Rust " + ISOLATED_RUST_VERSION)
        return compiler


def build_converter(tools_dir, rustc=None):
    """Build/cache only the std-only installation converter; never invoke Cargo."""
    source = complete_file(Path(__file__).resolve().with_name("spliceai2_converter.rs"))
    before = identity(source)
    source_sha = digest(source)
    require(identity(source) == before, "Installation converter source changed during hashing")
    tools_dir = Path(tools_dir).absolute()
    require(not any(parent.is_symlink() for parent in (tools_dir, *tools_dir.parents)),
            "Installation tool paths must not traverse symlinks")
    tools_dir.mkdir(parents=True, exist_ok=True)
    destination = tools_dir / source_sha
    executable = destination / "spliceai2-converter"
    with (tools_dir / ".build.lock").open("a") as lock:
        fcntl.flock(lock, fcntl.LOCK_EX)
        if destination.exists():
            receipt = destination / "build.ready.json"
            require(receipt.is_file() and executable.is_file() and os.access(executable, os.X_OK),
                    "Incomplete installation converter cache: " + str(destination))
            metadata = json.loads(receipt.read_text())
            require(metadata.get("schema_version") == 1 and metadata.get("source_sha256") == source_sha
                    and metadata.get("binary_sha256") == digest(executable),
                    "Installation converter cache identity mismatch: " + str(destination))
            return executable
        compiler = installation_compiler(tools_dir, rustc=rustc)
        compiler_inputs = {Path(compiler[key]): identity(Path(compiler[key])) for key in ("path", "executable")}
        with tempfile.TemporaryDirectory(prefix=".building-", dir=tools_dir) as temporary:
            stage = Path(temporary)
            binary = stage / executable.name
            command = [compiler["path"], "--edition=2021", "-O", "-C", "codegen-units=1", "--crate-name", "spliceai2_converter",
                       str(source), "-o", str(binary)]
            print("Building standalone SpliceAI2 installation converter", flush=True)
            subprocess.run(command, check=True)
            require(identity(source) == before, "Installation converter source changed during compilation")
            require(all(identity(path) == original for path, original in compiler_inputs.items()),
                    "Installation compiler changed during compilation")
            require(subprocess.check_output([str(binary), "--version"], text=True).strip()
                    == "spliceai2-converter 1 (installation-only; Rust std)",
                    "Unexpected installation converter identity")
            (stage / "build.ready.json").write_text(json.dumps({
                "schema_version": 1, "role": "installation_only", "source_path": str(source),
                "source_sha256": source_sha, "binary_sha256": digest(binary), "compiler": compiler,
                "command": command, "cargo_dependencies": [],
            }, indent=2) + "\n")
            stage.rename(destination)
    return executable


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


def verified_plugin(item, identities):
    path = verified_file(item, identities)
    sha256 = item["sha256"]
    require(sha256 in QUALIFIED_PLUGINS, "Unqualified official plugin revision")
    return {"path": str(path), "sha256": sha256, **QUALIFIED_PLUGINS[sha256]}


def text_stream(path):
    return gzip.open(path, "rt", encoding="utf-8", newline="") if path.suffix == ".gz" else path.open(encoding="utf-8", newline="")


def unversioned(gene_id):
    return re.sub(r"\.[0-9]+(_PAR_Y)?$", "", gene_id)


def vep_release_names(annotation, ensembl_release):
    """Gene names VEP reports, from the GENCODE GTF or GFF3 of VEP's Ensembl release."""
    gff3 = annotation.name.endswith((".gff3", ".gff3.gz"))
    require(gff3 or annotation.name.endswith((".gtf", ".gtf.gz")), "The VEP-release gene annotation must be a GTF or GFF3")
    names, description = {}, []
    with text_stream(annotation) as stream:
        for number, line in enumerate(stream, 1):
            if line.startswith("#"):
                description.append(line)
                continue
            fields = line.rstrip("\r\n").split("\t")
            require(len(fields) == 9, f"Invalid VEP-release annotation row {number}")
            if fields[2] != "gene":
                continue
            if gff3:
                attributes = {key: unquote(value) for key, value in (item.split("=", 1) for item in fields[8].split(";") if item)}
            else:
                attributes = dict(re.findall(r'(\w+) "([^"\n]*)";', fields[8]))
            gene_id, symbol = attributes.get("gene_id"), attributes.get("gene_name")
            require(gene_id and symbol, f"VEP-release annotation gene lacks ID/name at {number}")
            names.setdefault((unversioned(gene_id), fields[0]), set()).add((fields[6], symbol))
    description = "".join(description)
    require("GRCh38" in description and f"(Ensembl {ensembl_release})" in description,
            f"VEP gene names require the GRCh38 GENCODE annotation of Ensembl release {ensembl_release}")
    return names


def gene_context(regions, gtf, destination, release_names=None):
    """Resolve explicit Ensembl stable IDs, then require exact release/context."""
    genes, stable_ids, description = {}, {}, []
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
            match = re.fullmatch(r"(ENSG[0-9]{11})\.[0-9]+(_PAR_Y)?", gene_id)
            if match:
                stable_id = match[1] + (match[2] or "")
                stable_ids.setdefault(stable_id, []).append(gene_id)
    description = "".join(description)
    require("GRCh38" in description and "version 48" in description,
            "Gene mapping requires the original GRCh38 GENCODE version 48 GTF")
    rows, conventions, seen, resolved_seen, symbols = [], set(), set(), set(), {}
    stable_matches = 0
    with text_stream(regions) as stream:
        reader = csv.DictReader(stream, delimiter="\t")
        require(reader.fieldnames == REGION_COLUMNS, "Unexpected publisher scored-region schema")
        for row in reader:
            require(None not in row and all(value is not None for value in row.values()), "Malformed scored-region row")
            gene_id = row["gene_id"]
            require(gene_id not in seen, f"Duplicate scored-region gene: {gene_id}")
            seen.add(gene_id)
            resolved_id = gene_id
            if gene_id not in genes:
                candidates = stable_ids.get(gene_id, [])
                require(len(candidates) == 1,
                        f"No exact or unique GENCODE48 stable gene ID: {gene_id}; aliases/ambiguous versions are not allowed")
                resolved_id = candidates[0]
                stable_matches += 1
            require(resolved_id not in resolved_seen, f"Multiple source IDs resolve to one GENCODE gene: {resolved_id}")
            resolved_seen.add(resolved_id)
            chrom, start, end, strand, name, gene_type = genes[resolved_id]
            require(chrom == row["chrom"] and strand == row["strand"], f"Gene contig/strand mismatch: {gene_id}")
            require(gene_type == "protein_coding", f"Scored region is not a GENCODE48 protein-coding gene: {gene_id}")
            raw_start, raw_end = int(row["gene_start"]), int(row["gene_end"])
            require(raw_end == end and raw_start in (start, start - 1), f"Scored-region/GTF bounds differ: {gene_id}")
            conventions.add("1-based-closed" if raw_start == start else "0-based-half-open")
            require(chrom in CHROMS and strand in ("+", "-") and 0 < start <= end, f"Unsupported gene context: {gene_id}")
            require(not re.search(r"[\s|;,=]", name), f"Gene symbol cannot be represented safely in the plugin: {gene_id}")
            prior = symbols.setdefault((chrom, name), {"strand": strand, "ids": []})
            require(prior["strand"] == strand,
                    f"Opposite-strand genes share a symbol on {chrom}: {name}; plugin matching would be ambiguous")
            prior["ids"].append(resolved_id)
            rows.append((resolved_id, chrom, start, end, strand, name))
    require(rows and len(conventions) == 1, "Empty or mixed-coordinate scored-region table")
    genes_scored, release = len(rows), None
    if release_names is not None:
        # A second row with the same gene ID, interval and strand carries the VEP
        # release's name. A name another scored gene uses on that chromosome is
        # skipped, so no score moves to a different gene's transcripts.
        release = {"same_name": 0, "not_found": 0, "strand_differs": 0, "renamed": [], "skipped_name_conflicts": []}
        held, proposed = {(row[1], row[5]) for row in rows}, {}
        for row in rows:
            found = release_names.get((unversioned(row[0]), row[1]), set())
            if len(found) != 1:
                release["not_found"] += 1
                continue
            (strand, name), = found
            if strand != row[4]:
                release["strand_differs"] += 1
            elif name == row[5]:
                release["same_name"] += 1
            else:
                proposed.setdefault((row[1], name), []).append(row)
        for (chrom, name), owners in sorted(proposed.items()):
            if len(owners) > 1 or (chrom, name) in held or re.search(r"[\s|;,=]", name):
                release["skipped_name_conflicts"] += [[owner[0], owner[5], name] for owner in owners]
                continue
            rows.append(owners[0][:5] + (name,))
            release["renamed"].append([owners[0][0], owners[0][5], name])
    rows.sort(key=lambda row: (CHROMS.index(row[1]), row[2], row[3], row[0], row[5]))
    with destination.open("x", newline="") as stream:
        writer = csv.writer(stream, delimiter="\t", lineterminator="\n")
        writer.writerow(CONTEXT_COLUMNS)
        writer.writerows(rows)
    return {"genes": genes_scored, "context_rows": len(rows), "vep_release_names": release,
            "source_coordinates": conventions.pop(),
            "context_coordinates": "1-based-closed",
            "gene_id_matching": "unique-GENCODE48-stable-ID-with-exact-context" if stable_matches else "exact",
            "stable_id_matches": stable_matches,
            "versioned_gene_ids_retained": True,
            "same_strand_symbol_groups": sum(len(group["ids"]) > 1 for group in symbols.values()),
            "symbol_projection": "One unchanged prediction per distinct symbol per native row, under the GENCODE48 name and, when it differs, the VEP release's name of the same gene ID; all same-strand gene IDs/intervals retained in context. Opposite-strand same-symbol contexts rejected.",
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
                command = [str(converter), "project", "--genes", str(context), "--chrom", item["chrom"],
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
    require(isinstance(manifest.get("vep_gene_annotation"), dict),
            "A pinned GENCODE GTF/GFF3 of VEP's Ensembl release is required so renamed genes keep their scores")
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
    vep_annotation = verified_file(manifest["vep_gene_annotation"], identities)
    plugin = verified_plugin(manifest["vep_spliceai_plugin"], identities)
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
    if args.converter_bin:
        converter_name = shutil.which(str(args.converter_bin))
        require(converter_name, "The explicit --converter-bin installation tool is not executable")
        converter = complete_file(Path(converter_name).absolute())
        converter_info = {"role": "explicit_converter_override"}
    else:
        source = complete_file(driver.with_name("spliceai2_converter.rs"))
        identities[source] = identity(source)
        converter = build_converter(folder.parent / ".tools" / "spliceai2-converter")
        build_receipt = converter.parent / "build.ready.json"
        converter_info = {"role": "standalone_installation", "source_sha256": digest(source),
                          "build_receipt": str(build_receipt), "build_receipt_sha256": digest(build_receipt)}
    identities[converter] = identity(converter)
    converter_info.update(path=str(converter), sha256=digest(converter))
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
            # PriVA installs the plugin and the VEP cache from the same Ensembl release.
            mapping = gene_context(regions, gtf, context, vep_release_names(vep_annotation, plugin["release"]))
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
                "vep_gene_annotation": manifest["vep_gene_annotation"],
                "gene_context": {"path": str(folder / context.name), "sha256": digest(context), **mapping},
                "vep_spliceai_plugin": plugin,
                "converter": converter_info, "tools": tool_identities,
                "driver_sha256": driver_sha, "slurm_job_id": os.environ.get("SLURM_JOB_ID"),
                "thresholds": {"high_recall": 0.1, "balanced": 0.25, "high_precision": 0.5},
                "projection": "Unchanged top0 site delta/distance strings only, one allele/gene per VCF record. No cutoff filtering, score recomputation, masking, normalization or strand flipping. Original 84-column tables retain all top3/junction/reference/alternate/summary values.",
                "coverage_limits": "GENCODE48 protein-coding scored regions; population-observed indels; no published chrY indel shard; the plugin cannot score transcripts whose gene has no name in the VEP cache. Missing coverage is not zero score.",
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
    commands = parser.add_subparsers(dest="command", required=True)
    prepare_parser = commands.add_parser("prepare", help="Verify inputs and publish annotation-only VCFs")
    prepare_parser.add_argument("--manifest", type=Path, required=True)
    prepare_parser.add_argument("--assembly", required=True)
    prepare_parser.add_argument("--output-dir", type=Path, required=True)
    prepare_parser.add_argument("--receipt", type=Path, required=True)
    prepare_parser.add_argument("--converter-bin", help="Explicit standalone installation converter; default builds/caches the bundled source")
    prepare_parser.add_argument("--threads", type=int, default=1)
    prepare_parser.add_argument("--validation-subset", action="store_true", help="Explicit bounded test; never produces an installable full-cache receipt")
    build_parser = commands.add_parser("build-converter", help="Build or reuse the installation converter without processing score data")
    build_parser.add_argument("--tools-dir", type=Path, required=True)
    build_parser.add_argument("--rustc", help="Explicit compatible Rust compiler; default acquires an isolated compiler if needed")
    context_parser = commands.add_parser("gene-context", help="Resolve pinned scored regions against GENCODE48")
    context_parser.add_argument("--regions", type=Path, required=True)
    context_parser.add_argument("--gtf", type=Path, required=True)
    context_parser.add_argument("--output", type=Path, required=True)
    context_parser.add_argument("--vep-annotation", type=Path, help="GENCODE GTF/GFF3 of VEP's Ensembl release; adds its gene names")
    context_parser.add_argument("--ensembl-release", help="VEP's Ensembl release, required with --vep-annotation")
    argv = sys.argv[1:]
    # Preserve the previous prepare-only flags interface while exposing procedures.
    if argv and argv[0].startswith("--") and argv[0] != "--help":
        argv = ["prepare", *argv]
    args = parser.parse_args(argv)
    try:
        if args.command == "build-converter":
            print(build_converter(args.tools_dir, rustc=args.rustc))
        elif args.command == "gene-context":
            require((args.vep_annotation is None) == (args.ensembl_release is None),
                    "--vep-annotation and --ensembl-release go together")
            names = args.vep_annotation and vep_release_names(args.vep_annotation, args.ensembl_release)
            print(json.dumps(gene_context(args.regions, args.gtf, args.output, names)))
        else:
            prepare(args)
    except (ValueError, OSError, subprocess.SubprocessError, KeyError) as error:
        parser.exit(1, f"SpliceAI2 preparation failed: {error}\n")


if __name__ == "__main__":
    main()
