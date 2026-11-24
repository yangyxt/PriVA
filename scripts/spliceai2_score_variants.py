#!/usr/bin/env python3
"""Scores variants with the SpliceAI2 model on CPUs.

    spliceai2_score_variants.py --variants LIST.tsv --fasta GENOME.fa --model-dir FOLDER
        --output SCORES.tsv [--processes N] [--threads N]
    spliceai2_score_variants.py --version

LIST.tsv has the columns chrom, pos, ref, alt and strand (+ or -, the gene's strand). SCORES.tsv
gets one line per scored variant: the list's columns, then the strongest acceptor gain and loss
and donor gain and loss (DS_AG, DS_AL, DS_DG, DS_DL, two decimals) and their distances from the
variant (DP_AG, DP_AL, DP_DG, DP_DL). These are the values the published tables hold, and VEP's
SpliceAI plugin reads, under the same names.

The model, its 196,608-base window, the average over the weights in FOLDER/*/model.ckpt and the
extraction of the strongest changes are the publisher's spliceai2 package. Its own command needs a
GPU; this one runs the same prediction step on CPUs in float32. A variant whose window reaches
beyond a chromosome end is left out, because the package would score an empty sequence. A list
whose reference allele is not the genome's stops the run.

Rust lists the variants and writes the scores into the family table (rust/src/spliceai2_scores.rs).
"""
import argparse
import glob
import multiprocessing
import os
import sys
import time

WINDOW = 196608
PACKAGE_OUTPUT_LENGTH = 1024  # the package's default --out_2d_len
VARIANT_COLUMNS = ["chrom", "pos", "ref", "alt", "strand"]
EFFECTS = [("AG", "acceptor_gain"), ("AL", "acceptor_loss"), ("DG", "donor_gain"), ("DL", "donor_loss")]
STATE = {}


def start(model_dir, fasta, threads):
    """One process's model and genome."""
    os.environ["OMP_NUM_THREADS"] = str(threads)
    import pyfaidx
    import torch
    from spliceai2.model import SpliceAI2, VariantAnnotator

    torch.set_num_threads(threads)
    paths = sorted(glob.glob(os.path.join(model_dir, "*", "model.ckpt")))
    if not paths:
        raise SystemExit(f"{model_dir} holds no FOLDER/model.ckpt")
    models = [SpliceAI2.load_from_checkpoint(path, map_location="cpu").eval().float() for path in paths]
    STATE.update(
        torch=torch,
        model=VariantAnnotator(models, PACKAGE_OUTPUT_LENGTH).eval(),
        # The index beside the genome is used as it is; a shared genome is never written to.
        genome=pyfaidx.Fasta(fasta, rebuild=False),
        fasta=fasta,
    )


def window_inside(row):
    """Whether the package's window around the variant lies on the chromosome."""
    genome = STATE["genome"]
    if row["chrom"] not in genome:
        return False
    first = (row["pos"] - 1) - WINDOW // 2
    last = (row["pos"] - 1) + WINDOW // 2 + max(0, len(row["ref"]) - len(row["alt"]))
    return first >= 0 and last <= len(genome[row["chrom"]])


def score(item):
    """(index, values or None, seconds) of one variant."""
    import pandas
    from spliceai2.dataset import VariantDataset

    index, row = item
    began = time.time()
    if not window_inside(row):
        return index, None, time.time() - began
    torch = STATE["torch"]
    dataset = VariantDataset(pandas.DataFrame([dict(row, assembly="GRCh38")]), STATE["fasta"], WINDOW)
    dataset._fasta = STATE["genome"]
    with torch.inference_mode():
        result, columns = STATE["model"].predict_step(torch.utils.data.default_collate([dataset[0]]))
    value = dict(zip(columns, (float(number) for number in result[0])))
    scores = ["%.2f" % value[f"{name}_delta_score_0"] for _, name in EFFECTS]
    distances = [str(int(round(value[f"{name}_dist_0"]))) for _, name in EFFECTS]
    return index, scores + distances, time.time() - began


def read_variants(path):
    lines = open(path).read().splitlines()
    if not lines or lines[0].split("\t") != VARIANT_COLUMNS:
        raise SystemExit(f"{path} is not a list of variants ({', '.join(VARIANT_COLUMNS)})")
    rows = []
    for line in lines[1:]:
        fields = line.split("\t")
        if len(fields) != len(VARIANT_COLUMNS) or fields[4] not in ("+", "-"):
            raise SystemExit(f"{path}: a line is not a variant with its strand: {line!r}")
        rows.append(dict(zip(VARIANT_COLUMNS, fields), pos=int(fields[1])))
    return rows


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--version", action="store_true", help="print the spliceai2 package's version")
    parser.add_argument("--variants")
    parser.add_argument("--fasta")
    parser.add_argument("--model-dir")
    parser.add_argument("--output")
    parser.add_argument("--processes", type=int, default=1)
    parser.add_argument("--threads", type=int, default=1, help="CPUs per process")
    args = parser.parse_args()
    if args.version:
        from importlib.metadata import version

        print("spliceai2", version("spliceai2"))
        return
    for name in ("variants", "fasta", "model_dir", "output"):
        if not getattr(args, name):
            parser.error(f"--{name.replace('_', '-')} is required")
    rows = read_variants(args.variants)
    if not os.path.isfile(args.fasta + ".fai"):
        raise SystemExit(f"{args.fasta}.fai is missing; index the genome first")
    processes = max(1, min(args.processes, len(rows)))
    began = time.time()
    scored, seconds = {}, 0.0
    if rows:
        # Spawned, so that each process starts its own threads.
        context = multiprocessing.get_context("spawn")
        with context.Pool(processes, start, (args.model_dir, args.fasta, args.threads)) as pool:
            for index, values, spent in pool.imap_unordered(score, list(enumerate(rows)), chunksize=1):
                seconds += spent
                if values is not None:
                    scored[index] = values
    staged = args.output + ".partial"
    with open(staged, "w") as handle:
        handle.write("\t".join(VARIANT_COLUMNS + [f"DS_{k}" for k, _ in EFFECTS] + [f"DP_{k}" for k, _ in EFFECTS]) + "\n")
        for index, row in enumerate(rows):
            if index in scored:
                handle.write("\t".join([str(row[name]) for name in VARIANT_COLUMNS] + scored[index]) + "\n")
    os.replace(staged, args.output)
    print(
        f"SpliceAI2: {len(scored)} of {len(rows)} variants scored in {(time.time() - began) / 60:.1f} minutes"
        f" ({seconds / max(1, len(scored)):.1f} seconds each, {processes} processes of {args.threads} CPUs)",
        file=sys.stderr,
    )


if __name__ == "__main__":
    main()
