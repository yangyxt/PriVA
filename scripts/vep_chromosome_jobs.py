#!/usr/bin/env python3
"""Prepare and run independent official-VEP cache workers; never publish a cache.

User, 2026-10-23 CT: ten allocated CPUs per chromosome, up to twelve
chromosomes together, retaining completed work. In installed VEP116.3 the
Runner.pm spawn loop uses <=, so --fork N can create N+1 children in addition
to the parent. Eight forks therefore fit the proposed ten-CPU task.

The split is lossless, including duplicate, nonvariant and noncanonical rows.
VEP's own treatment of those rows is checked separately. A successful worker
means annotation finished, not output equivalence or permission to publish.
Use a Slurm array to call `worker` with independent ten-CPU allocations; the
caller must check live queue limits before choosing array concurrency. This
module neither submits jobs nor changes the input/configuration/canonical VCF.
"""
from __future__ import annotations

import argparse
from collections import Counter
import hashlib
import json
import os
from pathlib import Path
import re
import shutil
import subprocess

import yaml

from clinvar_gene_refresh import scientific_profile, resource_policy, profile_signature
from vep_chromosome_cache import (digest, identity, index_inventory, now, owned,
                                  open_vcf, read_header, require, row_invariant,
                                  same_identity, write_new, vep116_omission_policy,
                                  vep116_nonvariant_omitted)


KIND = "priva-vep-chromosome-worker-plan-v1"


def profile(config, source_root):
    return scientific_profile(config, source_root, resource_policy(source_root))


def pinned(path):
    return dict(identity=identity(path), sha256=digest(path))


def check_pin(value):
    path = Path(value["identity"]["path"])
    require(same_identity(value["identity"], identity(path)), "Pinned file identity changed: " + str(path))
    require(digest(path) == value["sha256"], "Pinned file checksum changed: " + str(path))
    require(same_identity(value["identity"], identity(path)), "Pinned file changed during checksum: " + str(path))
    return path


def read_plan(args):
    require(digest(args.manifest) == args.manifest_sha256, "Manifest checksum differs")
    plan = json.loads(args.manifest.read_text())
    require(plan.get("kind") == KIND, "Unsupported chromosome plan")
    require(plan["forks"] + 2 <= plan["cpus_per_chromosome"], "Worker CPU accounting is invalid")
    return plan


def split(args):
    require(args.forks > 0 and args.forks + 2 <= args.cpus_per_chromosome,
            "Allocate forks + two CPUs for this pinned VEP spawn loop")
    require(args.max_parallel > 0, "Positive concurrency required")
    source_root = args.source_root.absolute()
    config = yaml.safe_load(args.config.read_text())
    require(isinstance(config, dict) and config.get("cache_root"), "Configured cache_root is required")
    cache_root = Path(config["cache_root"])
    require(cache_root.is_absolute() and cache_root.is_dir(), "Existing absolute cache_root required")
    work = owned(args.work_dir, cache_root)
    require(work != cache_root and work.parent.is_dir() and not os.path.lexists(work),
            "Fresh work directory with an existing parent required")
    before_config = pinned(args.config)
    indexed = index_inventory(args.input_vcf, args.bcftools)
    worker = source_root / "scripts/annotation_vcf.sh"
    sources = {name: pinned(source_root / "scripts" / name)
               for name in ("annotation_vcf.sh", "vep_chromosome_cache.py", "vep_chromosome_jobs.py",
                            "clinvar_gene_refresh.py", "argparse.bash", "anno_vep_args")}
    science = profile(config, source_root)
    work.mkdir(mode=0o700)
    frozen_config = work / "config.yaml"
    with args.config.open("rb") as src, frozen_config.open("xb") as out:
        shutil.copyfileobj(src, out)
    require(digest(frozen_config) == before_config["sha256"], "Configuration copy differs")
    check_pin(before_config)
    frozen_config.chmod(0o400)
    counts, shards = Counter(), []
    header = None
    current = None
    output = None
    previous_rank, previous_position = -1, 0
    order = indexed["contig_order"]
    ranks = {contig: i for i, contig in enumerate(order)}
    raw_body_hash = hashlib.sha256()
    ordinal = 0
    try:
        with open_vcf(args.input_vcf) as src:
            header = read_header(src)
            require(hashlib.sha256(header).hexdigest() == indexed["header_sha256"], "Source header changed")
            for line in src:
                invariant = row_invariant(line)
                contig = invariant[0].decode()
                require(contig in ranks, "Source contig absent from its index")
                rank, position = ranks[contig], int(invariant[1])
                require(rank >= previous_rank and (rank != previous_rank or position >= previous_position),
                        "Source contig/position order regressed")
                previous_rank, previous_position = rank, position
                if contig != current:
                    if output is not None:
                        output.close()
                        shards[-1].update(input=pinned(shards[-1]["input_path"]),
                                          body_sha256=shard_hash.hexdigest())
                        require(shards[-1]["input"]["sha256"] == shard_full_hash.hexdigest(),
                                "Chromosome copy bytes differ from the source stream")
                    folder = work / f"contig-{rank:03d}"
                    folder.mkdir(mode=0o700)
                    path = folder / "input.vcf"
                    output = path.open("xb")
                    output.write(header)
                    shards.append(dict(index=rank, contig=contig, input_path=str(path),
                                       source_ordinal_start=ordinal + 1, folder=str(folder)))
                    shard_hash = hashlib.sha256()
                    shard_full_hash = hashlib.sha256(header)
                    current = contig
                output.write(line)
                raw_body_hash.update(line)
                shard_hash.update(line)
                shard_full_hash.update(line)
                counts[contig] += 1
                ordinal += 1
            if output is not None:
                output.close()
                output = None
                shards[-1].update(input=pinned(shards[-1]["input_path"]), body_sha256=shard_hash.hexdigest())
                require(shards[-1]["input"]["sha256"] == shard_full_hash.hexdigest(),
                        "Chromosome copy bytes differ from the source stream")
        require(dict(counts) == indexed["records_by_contig"], "Split source/index record counts differ")
        require(same_identity(indexed["input"], identity(args.input_vcf)), "Source changed during split")
        require(profile_signature(profile(config, source_root)) == profile_signature(science),
                "Scientific resources changed during split")
        for value in sources.values():
            check_pin(value)
        for shard in shards:
            shard["source_records"] = counts[shard["contig"]]
            shard["source_ordinal_end"] = shard["source_ordinal_start"] + shard["source_records"] - 1
            Path(shard["input_path"]).chmod(0o400)
            shard["input"]["identity"] = identity(shard["input_path"])
        plan = dict(kind=KIND, created_utc=now(), work_dir=str(work), cache_root=str(cache_root),
                    source_root=str(source_root), worker=str(worker), sources=sources,
                    original_config=before_config, worker_config=pinned(frozen_config),
                    scientific_profile=science,
                    input_inventory=indexed, source_body_sha256=raw_body_hash.hexdigest(),
                    source_records=ordinal, shards=shards, forks=args.forks,
                    cpus_per_chromosome=args.cpus_per_chromosome, max_parallel=args.max_parallel,
                    upstream_record_filtering_applied_during_split=False,
                    canonical_publication=False, installation_accepted=False, scientific_qualification=False)
        manifest = work / "chromosomes.json"
        write_new(manifest, plan)
        manifest.chmod(0o400)
        print(json.dumps(dict(manifest=str(manifest), sha256=digest(manifest), chromosomes=len(shards),
                              source_records=ordinal, annotated_completed_chromosomes=[])))
    finally:
        if output is not None:
            output.close()


def worker(args):
    plan = read_plan(args)
    require(0 <= args.contig_index < len(plan["shards"]), "Chromosome index out of range")
    shard = plan["shards"][args.contig_index]
    require(shard["index"] == args.contig_index, "Chromosome index/order mismatch")
    for value in plan["sources"].values():
        check_pin(value)
    config = check_pin(plan["worker_config"])
    source = check_pin(shard["input"])
    current_profile = profile(yaml.safe_load(config.read_text()), Path(plan["source_root"]))
    require(profile_signature(current_profile) == profile_signature(plan["scientific_profile"]),
            "Worker scientific resources differ from the split plan")
    cpus = os.environ.get("SLURM_CPUS_PER_TASK")
    require(cpus is not None and cpus.isdigit() and int(cpus) >= plan["cpus_per_chromosome"],
            "Worker requires its declared Slurm CPU allocation")
    root = owned(shard["folder"], plan["cache_root"])
    output, stats, log = root / "output.vcf", root / "stats.html", root / "worker.log"
    started, receipt = root / "started.json", root / "worker-result.json"
    require(not any(os.path.lexists(p) for p in (output, stats, log, started, receipt)),
            "Worker attempt already exists; preserve it and prepare a fresh retry")
    command = ["bash", plan["worker"], "anno_VEP_cache_worker", "--config", str(config),
               "--input-vcf", str(source), "--output-vcf", str(output), "--stats-file", str(stats),
               "--forks", str(plan["forks"])]
    context = dict(kind="priva-vep-chromosome-worker-result-v1", contig=shard["contig"],
                   manifest=str(args.manifest.absolute()), manifest_sha256=args.manifest_sha256,
                   started_utc=now(), argv=command,
                   slurm={key: os.environ.get(key) for key in
                          ("SLURM_JOB_ID", "SLURM_ARRAY_JOB_ID", "SLURM_ARRAY_TASK_ID",
                           "SLURM_STEP_ID", "SLURM_CPUS_PER_TASK", "SLURMD_NODENAME")},
                   output_equivalence_validated=False, canonical_publication=False)
    write_new(started, context)
    with log.open("xb") as handle:
        result = subprocess.run(command, stdin=subprocess.DEVNULL, stdout=handle, stderr=subprocess.STDOUT)
    validation_error = None
    try:
        after_profile = profile(yaml.safe_load(config.read_text()), Path(plan["source_root"]))
        unchanged = profile_signature(after_profile) == profile_signature(current_profile)
    except Exception as error:
        # A closed child with a failed postcondition is failed evidence, not
        # a supposedly still-running worker. Preserve its actual exit code.
        unchanged = False
        validation_error = type(error).__name__ + ": " + str(error)
    context.update(completed_utc=now(), returncode=result.returncode,
                   scientific_resources_unchanged=unchanged,
                   postcondition_error=validation_error,
                   state="annotation_finished_pending_validation" if result.returncode == 0 and unchanged else "failed",
                   log=pinned(log))
    for label, path in (("output", output), ("stats", stats)):
        if path.is_file() and not path.is_symlink():
            # The annotation body is hashed by the later streaming validator,
            # avoiding an extra full VCF hash here. No ready claim is made.
            context[label] = identity(path)
    write_new(receipt, context)
    print(json.dumps(dict(contig=shard["contig"], state=context["state"], receipt=str(receipt))))
    if result.returncode:
        raise SystemExit(result.returncode if result.returncode > 0 else 1)
    require(unchanged, "Worker scientific resources changed while annotating")


def status(args):
    plan = read_plan(args)
    rows = []
    for shard in plan["shards"]:
        root = Path(shard["folder"])
        path = root / "worker-result.json"
        state = "started_not_closed" if (root / "started.json").exists() else "pending"
        if path.exists():
            result = json.loads(path.read_text())
            require(result.get("manifest_sha256") == args.manifest_sha256
                    and result.get("contig") == shard["contig"], "Worker receipt belongs to another plan")
            state = result["state"]
        proof = root / "source-occurrences.validated.json"
        if proof.exists():
            validated = json.loads(proof.read_text())
            require(validated.get("manifest_sha256") == args.manifest_sha256
                    and validated.get("contig") == shard["contig"]
                    and validated.get("source_occurrences_validated") is True,
                    "Validation proof belongs to another chromosome plan")
            require(same_identity(validated["output"]["identity"], identity(root / "output.vcf")),
                    "Validated chromosome output changed")
            state = "source_occurrences_validated_annotations_pending_equivalence"
        rows.append(dict(contig=shard["contig"], state=state, source_records=shard["source_records"]))
    print(json.dumps(dict(chromosomes=rows, canonical_publication=False, installation_accepted=False)))


def validate(args):
    """Account for every source occurrence without claiming annotation parity.

    This is distinct from the cross-profile real VEP comparison. The latter
    must succeed before a changed fork/chromosome execution route is adopted.
    """
    plan = read_plan(args)
    require(0 <= args.contig_index < len(plan["shards"]), "Chromosome index out of range")
    shard = plan["shards"][args.contig_index]
    for value in plan["sources"].values():
        check_pin(value)
    config = check_pin(plan["worker_config"])
    science = profile(yaml.safe_load(config.read_text()), Path(plan["source_root"]))
    require(profile_signature(science) == profile_signature(plan["scientific_profile"]),
            "Validation scientific resources differ from the split plan")
    source = check_pin(shard["input"])
    folder = owned(shard["folder"], plan["cache_root"])
    result_path = folder / "worker-result.json"
    before_result = pinned(result_path)
    result = json.loads(result_path.read_text())
    require(result.get("kind") == "priva-vep-chromosome-worker-result-v1"
            and result.get("manifest_sha256") == args.manifest_sha256
            and result.get("contig") == shard["contig"]
            and result.get("state") == "annotation_finished_pending_validation"
            and result.get("returncode") == 0 and result.get("scientific_resources_unchanged") is True,
            "A matching successful closed VEP child is required")
    output = folder / "output.vcf"
    before_output = identity(output)
    require(same_identity(before_output, result["output"]), "Worker output changed since its child closed")
    omissions = folder / "validation-omitted-source-occurrences.jsonl"
    proof = folder / "source-occurrences.validated.json"
    require(not os.path.lexists(omissions) and not os.path.lexists(proof), "Validation attempt already exists")
    kept, skipped, total = 0, 0, 0
    ledger_hash, body_hash = hashlib.sha256(), hashlib.sha256()
    with open_vcf(source) as src, output.open("rb") as annotated, omissions.open("xb") as ledger:
        source_header, output_header = read_header(src), read_header(annotated)
        require(hashlib.sha256(source_header).hexdigest() == plan["input_inventory"]["header_sha256"],
                "Chromosome source header differs from the split plan")
        require(source_header.splitlines()[-1] == output_header.splitlines()[-1], "VCF column/sample header differs")
        output_lines = Counter(output_header.splitlines(keepends=True))
        for line, count in Counter(source_header.splitlines(keepends=True)).items():
            if line.startswith((b"##fileformat=", b"##INFO=", b"##contig=")):
                require(output_lines[line] == count, "Original scientific header changed")
        require(len(re.findall(rb"^##INFO=<ID=CSQ,", output_header, re.M)) == 1,
                "Missing/duplicate CSQ schema")
        policy = vep116_omission_policy(output_header, args.vep_source_root)
        output_hash = hashlib.sha256(output_header)
        source_hash = hashlib.sha256(source_header)
        for line in src:
            total += 1
            source_hash.update(line)
            invariant = row_invariant(line)
            require(invariant[0].decode() == shard["contig"], "Unexpected chromosome in split input")
            if policy["nonvariant_omissions_authorized"] and vep116_nonvariant_omitted(
                    invariant, allow_non_variant=policy["allow_non_variant"]):
                value = dict(source_ordinal=shard["source_ordinal_start"] + total - 1,
                             source_contig_ordinal=total, next_output_contig_ordinal=kept + 1,
                             source_row_sha256=hashlib.sha256(line).hexdigest(), source_row_bytes=len(line),
                             reason="vep116_first_alt_dot_allow_non_variant_false_nonstructural")
                encoded = (json.dumps(value, sort_keys=True, allow_nan=False) + "\n").encode()
                ledger.write(encoded)
                ledger_hash.update(encoded)
                skipped += 1
                continue
            row = annotated.readline()
            require(row_invariant(row, True) == invariant, "Missing/extra/reordered/changed output occurrence")
            output_hash.update(row)
            body_hash.update(row)
            kept += 1
        require(not annotated.read(1), "Extra output occurrence after all accepted source rows")
    require(total == shard["source_records"] and kept + skipped == total, "Source occurrence accounting differs")
    require(source_hash.hexdigest() == shard["input"]["sha256"], "Split source bytes changed during validation")
    require(same_identity(shard["input"]["identity"], identity(source)), "Split source identity changed")
    require(same_identity(before_output, identity(output)), "Worker output changed during validation")
    check_pin(before_result)
    ledger_pin = pinned(omissions)
    require(ledger_pin["sha256"] == ledger_hash.hexdigest(), "Omission ledger copy differs")
    require(profile_signature(profile(yaml.safe_load(config.read_text()), Path(plan["source_root"])))
            == profile_signature(science), "Scientific resources changed during validation")
    value = dict(kind="priva-vep-chromosome-source-occurrences-v1", created_utc=now(),
                 manifest=str(args.manifest.absolute()), manifest_sha256=args.manifest_sha256,
                 contig=shard["contig"], worker_result=before_result,
                 output=dict(identity=before_output, sha256=output_hash.hexdigest()),
                 output_body_sha256=body_hash.hexdigest(), output_header_sha256=hashlib.sha256(output_header).hexdigest(),
                 source_records=total, emitted_records=kept, upstream_omitted_records=skipped,
                 omitted_occurrences=ledger_pin, parser_omission_policy=policy,
                 source_ordinal_start=shard["source_ordinal_start"], source_ordinal_end=shard["source_ordinal_end"],
                 source_occurrences_validated=True, annotation_equivalence_validated=False,
                 canonical_publication=False, installation_accepted=False, scientific_qualification=False)
    write_new(proof, value)
    print(json.dumps(dict(proof=str(proof), sha256=digest(proof), contig=shard["contig"],
                          source_records=total, emitted_records=kept, upstream_omitted_records=skipped,
                          annotation_equivalence_validated=False)))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest="operation", required=True)
    p = sub.add_parser("split", help="Losslessly freeze config and partition the indexed source once")
    for name in ("config", "source-root", "input-vcf", "work-dir"):
        p.add_argument("--" + name, required=True, type=Path)
    p.add_argument("--bcftools", default="bcftools")
    p.add_argument("--forks", default=8, type=int)
    p.add_argument("--cpus-per-chromosome", default=10, type=int)
    p.add_argument("--max-parallel", default=12, type=int)
    p.set_defaults(func=split)
    for name, function in (("worker", worker), ("status", status), ("validate", validate)):
        p = sub.add_parser(name)
        p.add_argument("--manifest", required=True, type=Path)
        p.add_argument("--manifest-sha256", required=True)
        if name in ("worker", "validate"):
            p.add_argument("--contig-index", required=True, type=int)
        if name == "validate":
            p.add_argument("--vep-source-root", type=Path)
        p.set_defaults(func=function)
    args = parser.parse_args()
    args.func(args)


if __name__ == "__main__":
    main()
