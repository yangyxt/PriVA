#!/usr/bin/env python3
"""
Missense-Constrained Segment Extraction via Anchor-Guided Cross-Validated KDE.

Replaces the old min/max KDE approach with a principled workflow:
  1. Parse AlphaMissense VCF → per-residue composite score  x_i = f_i · μ_i · 100
     where f_i = fraction of substitutions with AM > 0.564 ("likely pathogenic"),
           μ_i = mean AM score across all substitutions at that position.
  2. Load positive anchors  (RMC constrained + HCSeeker hotspots)
     and   negative anchors (RMC unconstrained + HCSeeker coldspots).
  3. 2D cross-validation over (bandwidth h, threshold t) to jointly select h* and t*
     that maximise   Q(h,t) = Precision(predicted, positive) − η · Spillover(predicted, negative).
  4. Threshold candidates span from 1%-percentile of pos anchors to 99%-percentile of neg anchors.
  5. Extract contiguous constrained segments for every protein using (h*, t*).

Reference: missense_constrained_segments_plan.md
"""

import os
# Limit BLAS threads before importing numpy (must precede numpy import)
os.environ["OPENBLAS_NUM_THREADS"] = "2"
os.environ["MKL_NUM_THREADS"] = "2"
os.environ["NUMEXPR_NUM_THREADS"] = "2"
os.environ["OMP_NUM_THREADS"] = "2"

import re
import warnings
import sys
import time
import pysam
import pickle
import logging
import argparse
import numpy as np
import pandas as pd
import multiprocessing as mp
from pathlib import Path
from collections import defaultdict
from functools import partial
from sklearn.neighbors import KernelDensity
from sklearn.discriminant_analysis import LinearDiscriminantAnalysis
from sklearn.metrics import roc_auc_score
from scipy.stats import false_discovery_control, percentileofscore
from concurrent.futures import ProcessPoolExecutor, as_completed

# ──────────────────────────────────────── constants ────────────────────────────
AM_PATHOGENIC_THRESHOLD = 0.564  # AlphaMissense "likely pathogenic" cutoff (fixed)
COMPOSITE_SCALE = 100            # scaling factor for composite score

# Three-letter → one-letter amino acid code (for parsing RMC start/end_aa_raw)
AA3_TO_1 = {
    'Ala': 'A', 'Arg': 'R', 'Asn': 'N', 'Asp': 'D', 'Cys': 'C',
    'Gln': 'Q', 'Glu': 'E', 'Gly': 'G', 'His': 'H', 'Ile': 'I',
    'Leu': 'L', 'Lys': 'K', 'Met': 'M', 'Phe': 'F', 'Pro': 'P',
    'Ser': 'S', 'Thr': 'T', 'Trp': 'W', 'Tyr': 'Y', 'Val': 'V',
    'Sec': 'U', 'Pyl': 'O', 'Ter': '*',
}


def setup_logging(log_file: str | None = None) -> logging.Logger:
    """Configure root logger with both console and optional file output."""
    logger = logging.getLogger("mcs")  # missense_constrained_segments
    logger.setLevel(logging.DEBUG)
    fmt = logging.Formatter(
        "%(levelname)s:%(asctime)s:%(funcName)s:%(lineno)s: %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
    )
    # Console handler (INFO level)
    ch = logging.StreamHandler()
    ch.setLevel(logging.INFO)
    ch.setFormatter(fmt)
    logger.addHandler(ch)
    # File handler (DEBUG level) — captures everything
    if log_file:
        fh = logging.FileHandler(log_file, mode="w")
        fh.setLevel(logging.DEBUG)
        fh.setFormatter(fmt)
        logger.addHandler(fh)
    return logger


# ═══════════════════════════════════════════════════════════════════════════════
#  Section 1: VCF parsing → per-residue composite scores
# ═══════════════════════════════════════════════════════════════════════════════

def _parse_chromosome(args_tuple):
    """Parse a single chromosome from the VCF.  Designed for multiprocessing.

    Returns list of (transcript_id, [(aa_pos, ref_aa, am_score), ...]).
    """
    vcf_path, chromosome, gene_field, transcript_field, consq_field = args_tuple
    chrom_variants = defaultdict(list)
    count = 0
    try:
        vcf = pysam.VariantFile(vcf_path)
        for record in vcf.fetch(chromosome):
            try:
                am_score = float(record.info["AM_PATHOGENICITY"])
                target_tx = record.info.get("TRANSCRIPT").split(".")[0]
                pvar = record.info.get("PVAR")
                csqs = record.info.get("CSQ", [])
                if not csqs:
                    continue
                for csq in csqs:
                    fields = csq.split("|")
                    if fields[consq_field] != "missense_variant":
                        continue
                    if fields[transcript_field] != target_tx:
                        continue
                    if not fields[gene_field] or not pvar:
                        continue
                    m = re.match(r"^([A-Z]+)(\d+)([A-Z]+)$", pvar)
                    if not m:
                        continue
                    ref_aa, aa_pos = m.group(1), int(m.group(2))
                    chrom_variants[target_tx].append((aa_pos, ref_aa, am_score))
                    count += 1
            except (KeyError, ValueError):
                continue
        vcf.close()
    except Exception as e:
        logging.getLogger("mcs").error(f"Error on chr {chromosome}: {e}")
    return list(chrom_variants.items())


def parse_vcf(vcf_path: str, n_processes: int = 4) -> dict[str, list[tuple]]:
    """Parse AlphaMissense VEP-annotated VCF into per-transcript variant lists.

    Returns:
        dict  { transcript_id: [(aa_pos, ref_aa, am_score), ...] }
    """
    logger = logging.getLogger("mcs")
    logger.info(f"Parsing VCF: {vcf_path}")

    vcf = pysam.VariantFile(vcf_path)
    csq_fields = vcf.header.info["CSQ"].description.split("Format: ")[1].split("|")
    gene_idx = csq_fields.index("Gene")
    tx_idx = csq_fields.index("Feature")
    csq_idx = csq_fields.index("Consequence")
    chromosomes = list(vcf.header.contigs)
    vcf.close()

    logger.info(f"  {len(chromosomes)} chromosomes, {n_processes} workers")
    task_args = [(vcf_path, chrom, gene_idx, tx_idx, csq_idx) for chrom in chromosomes]

    with mp.Pool(n_processes) as pool:
        chrom_results = pool.map(_parse_chromosome, task_args)

    # Merge per-chromosome results
    transcript_variants: dict[str, list] = defaultdict(list)
    total = 0
    for items in chrom_results:
        for tx, variants in items:
            transcript_variants[tx].extend(variants)
            total += len(variants)

    logger.info(f"  Parsed {total:,} variant entries across {len(transcript_variants):,} transcripts")
    return dict(transcript_variants)


def compute_composite_scores(
    transcript_variants: dict[str, list[tuple]],
) -> dict[str, tuple[np.ndarray, np.ndarray, dict[int, str]]]:
    """Convert per-variant AM scores to per-residue composite scores.

    For each position i:
        f_i = fraction of substitutions with AM > 0.564
        μ_i = mean AM score
        x_i = f_i · μ_i · 100

    Returns:
        dict { transcript_id: (positions_array, composite_scores_array, {pos: ref_aa}) }
    """
    logger = logging.getLogger("mcs")
    result = {}
    for tx, variants in transcript_variants.items():
        if len(variants) < 5:  # too few data points for meaningful KDE
            continue
        # Group scores by position
        pos_scores: dict[int, list[float]] = defaultdict(list)
        pos_ref_aa: dict[int, str] = {}
        for pos, ref_aa, score in variants:
            pos_scores[pos].append(score)
            pos_ref_aa[pos] = ref_aa

        positions = np.array(sorted(pos_scores.keys()))
        n_pos = len(positions)
        composite = np.empty(n_pos, dtype=np.float64)

        for idx, p in enumerate(positions):
            scores_arr = np.array(pos_scores[p])
            f_i = np.mean(scores_arr > AM_PATHOGENIC_THRESHOLD)  # fraction pathogenic
            mu_i = scores_arr.mean()                              # mean AM
            composite[idx] = f_i * mu_i * COMPOSITE_SCALE

        result[tx] = (positions, composite, pos_ref_aa)

    logger.info(f"  Composite scores computed for {len(result):,} transcripts")
    return result


def compute_extended_features(
    transcript_variants: dict[str, list[tuple]],
) -> dict[str, tuple[np.ndarray, dict[str, np.ndarray], dict[int, str]]]:
    """Compute multiple per-residue features for each transcript.

    Features:
        f        fraction of substitutions with AM > 0.564
        mu       mean AM score
        min      min AM score
        max      max AM score
        median   median AM score
        std      std dev of AM scores
        n        number of substitutions at position
        f_mu     f × mu  (current composite / 100)
        f_max    f × max
        f_sq_mu  f² × mu

    Returns:
        dict { tx: (positions_array, {name: feature_array}, {pos: ref_aa}) }
    """
    logger = logging.getLogger("mcs")
    result = {}
    for tx, variants in transcript_variants.items():
        if len(variants) < 5:
            continue
        pos_scores: dict[int, list[float]] = defaultdict(list)
        pos_ref_aa: dict[int, str] = {}
        for pos, ref_aa, score in variants:
            pos_scores[pos].append(score)
            pos_ref_aa[pos] = ref_aa

        positions = np.array(sorted(pos_scores.keys()))
        n_pos = len(positions)

        features: dict[str, np.ndarray] = {
            k: np.empty(n_pos, dtype=np.float64)
            for k in ['f', 'mu', 'min', 'max', 'median', 'std', 'n',
                       'f_mu', 'f_max', 'f_sq_mu']
        }

        for idx, p in enumerate(positions):
            scores = np.array(pos_scores[p])
            f_i = np.mean(scores > AM_PATHOGENIC_THRESHOLD)
            mu_i = scores.mean()
            features['f'][idx] = f_i
            features['mu'][idx] = mu_i
            features['min'][idx] = scores.min()
            features['max'][idx] = scores.max()
            features['median'][idx] = np.median(scores)
            features['std'][idx] = scores.std() if len(scores) > 1 else 0.0
            features['n'][idx] = len(scores)
            features['f_mu'][idx] = f_i * mu_i
            features['f_max'][idx] = f_i * scores.max()
            features['f_sq_mu'][idx] = f_i ** 2 * mu_i

        result[tx] = (positions, features, pos_ref_aa)

    logger.info(f"  Extended features computed for {len(result):,} transcripts")
    return result


def audit_feature_separation(
    extended_data: dict[str, tuple[np.ndarray, dict[str, np.ndarray], dict]],
    pos_anchors: dict[str, set[int]],
    neg_anchors: dict[str, set[int]],
) -> None:
    """Measure AUC-ROC for each per-residue feature at anchor positions.

    Also fits Fisher LDA on all features to find the optimal linear combination
    and reports its AUC.
    """
    logger = logging.getLogger("mcs")

    # Collect feature values at anchor residues
    feature_names = None
    pos_rows: list[dict[str, float]] = []
    neg_rows: list[dict[str, float]] = []

    for tx, (positions, features, _) in extended_data.items():
        if feature_names is None:
            feature_names = list(features.keys())

        # Build a fast lookup: position → index
        pos_to_idx = {int(p): i for i, p in enumerate(positions)}

        for anchor_dict, target in [(pos_anchors, pos_rows), (neg_anchors, neg_rows)]:
            if tx not in anchor_dict:
                continue
            for apos in anchor_dict[tx]:
                idx = pos_to_idx.get(apos)
                if idx is None:
                    continue  # no variant data at this anchor position
                target.append({fn: features[fn][idx] for fn in feature_names})

    n_pos = len(pos_rows)
    n_neg = len(neg_rows)
    logger.info(f"Feature audit: {n_pos:,} positive anchor residues, "
                f"{n_neg:,} negative anchor residues with variant data")

    if n_pos < 10 or n_neg < 10:
        logger.error("Too few anchor residues for meaningful audit")
        return

    # Convert to arrays
    pos_arrays = {fn: np.array([r[fn] for r in pos_rows]) for fn in feature_names}
    neg_arrays = {fn: np.array([r[fn] for r in neg_rows]) for fn in feature_names}

    # ── Individual feature AUCs ──
    y_true = np.concatenate([np.ones(n_pos), np.zeros(n_neg)])
    results = []

    for fn in feature_names:
        p = pos_arrays[fn]
        n = neg_arrays[fn]
        y_score = np.concatenate([p, n])
        auc = roc_auc_score(y_true, y_score)
        # Flip if AUC < 0.5 (feature separates in opposite direction)
        if auc < 0.5:
            auc = 1.0 - auc
        pooled_std = np.sqrt((p.var() + n.var()) / 2)
        cohen_d = (p.mean() - n.mean()) / pooled_std if pooled_std > 0 else 0.0
        results.append((fn, auc, cohen_d, p.mean(), n.mean(), p.std(), n.std()))

    results.sort(key=lambda x: x[1], reverse=True)

    # ── Fisher LDA on all features ──
    X_pos = np.column_stack([pos_arrays[fn] for fn in feature_names])
    X_neg = np.column_stack([neg_arrays[fn] for fn in feature_names])
    X = np.vstack([X_pos, X_neg])
    y = np.concatenate([np.ones(len(X_pos)), np.zeros(len(X_neg))])

    lda = LinearDiscriminantAnalysis(n_components=1)
    lda.fit(X, y)
    lda_scores = lda.transform(X).ravel()
    lda_auc = roc_auc_score(y, lda_scores)

    # ── Log results ──
    logger.info("=" * 100)
    logger.info("FEATURE SEPARATION AUDIT  (raw per-residue, no KDE smoothing)")
    logger.info("=" * 100)
    logger.info(f"{'Feature':>12s} │ {'AUC':>7s} │ {'Cohen_d':>8s} │ "
                f"{'mu_pos':>9s} │ {'mu_neg':>9s} │ {'sd_pos':>9s} │ {'sd_neg':>9s}")
    logger.info("─" * 100)
    best_auc = results[0][1]
    for fn, auc, d, mp, mn, sp, sn in results:
        marker = " ◄ best" if auc == best_auc else ""
        logger.info(f"{fn:>12s} │ {auc:7.4f} │ {d:+8.4f} │ "
                    f"{mp:9.5f} │ {mn:9.5f} │ {sp:9.5f} │ {sn:9.5f}{marker}")
    logger.info("─" * 100)
    logger.info(f"{'LDA_optimal':>12s} │ {lda_auc:7.4f} │")
    logger.info("─" * 100)

    # LDA coefficients
    coefs = dict(zip(feature_names, lda.coef_[0]))
    coef_str = ", ".join(f"{fn}={w:+.6f}" for fn, w in
                         sorted(coefs.items(), key=lambda x: abs(x[1]), reverse=True))
    logger.info(f"LDA weights (sorted by |w|): {coef_str}")
    logger.info(f"LDA intercept: {lda.intercept_[0]:.6f}")
    logger.info("=" * 100)


# ═══════════════════════════════════════════════════════════════════════════════
#  Section 2: Anchor loading (RMC + HCSeeker)
# ═══════════════════════════════════════════════════════════════════════════════

def _parse_rmc_aa_pos(raw: str) -> int:
    """Parse RMC amino acid field like 'Arg97' → 97."""
    m = re.match(r"^([A-Za-z]+)(\d+)$", raw)
    if not m:
        raise ValueError(f"Cannot parse RMC AA position: {raw}")
    return int(m.group(2))


def load_rmc_anchors(
    rmc_path: str,
    oe_constrained: float = 0.4,
    oe_unconstrained: float = 0.8,
    p_bh_threshold: float = 0.05,
) -> tuple[dict[str, set[int]], dict[str, set[int]]]:
    """Load RMC regions and classify into positive / negative anchors.

    Applies Benjamini-Hochberg FDR correction to the raw p-values first.

    Constrained (positive):   oe < oe_constrained  AND  p_bh < p_bh_threshold
    Unconstrained (negative): oe > oe_unconstrained  OR  p_bh > p_bh_threshold
    Indeterminate:            everything else → excluded

    Returns:
        (pos_anchors, neg_anchors)  each is  dict { ENST_id: set of residue positions }
    """
    logger = logging.getLogger("mcs")
    logger.info(f"Loading RMC: {rmc_path}")

    df = pd.read_csv(rmc_path, sep="\t")
    logger.info(f"  {len(df)} regions, {df['transcript'].nunique()} transcripts")

    # Parse amino acid range endpoints
    df["start_pos"] = df["start_aa_raw"].apply(_parse_rmc_aa_pos)
    df["end_pos"] = df["end_aa_raw"].apply(_parse_rmc_aa_pos)

    # BH correction on raw p-values
    raw_p = df["p"].values.astype(np.float64)
    # scipy false_discovery_control returns adjusted p-values (Benjamini-Hochberg)
    df["p_bh"] = false_discovery_control(raw_p, method="bh")

    # Classify regions
    constrained_mask = (df["oe"] < oe_constrained) & (df["p_bh"] < p_bh_threshold)
    unconstrained_mask = (df["oe"] > oe_unconstrained) | (df["p_bh"] > p_bh_threshold)

    logger.info(
        f"  Constrained: {constrained_mask.sum()}, "
        f"Unconstrained: {unconstrained_mask.sum()}, "
        f"Indeterminate: {(~constrained_mask & ~unconstrained_mask).sum()}"
    )

    pos_anchors: dict[str, set[int]] = defaultdict(set)
    neg_anchors: dict[str, set[int]] = defaultdict(set)

    for mask, target in [(constrained_mask, pos_anchors), (unconstrained_mask, neg_anchors)]:
        for _, row in df[mask].iterrows():
            tx = row["transcript"]
            residues = set(range(row["start_pos"], row["end_pos"] + 1))
            target[tx].update(residues)

    logger.info(
        f"  Positive anchor transcripts: {len(pos_anchors)}, "
        f"Negative anchor transcripts: {len(neg_anchors)}"
    )
    return dict(pos_anchors), dict(neg_anchors)


def _load_nm_to_enst_map(mapping_path: str) -> dict[str, list[str]]:
    """Load RefSeq NM → Ensembl ENST mapping table.

    Returns dict { NM_id_no_version: [ENST_id, ...] }
    """
    df = pd.read_csv(mapping_path, sep="\t")
    # Column names: Gene stable ID, Transcript stable ID, HGNC symbol,
    #               UniProtKB/Swiss-Prot ID, RefSeq mRNA ID
    col_enst = df.columns[1]  # Transcript stable ID
    col_nm = df.columns[4]    # RefSeq mRNA ID
    df = df.dropna(subset=[col_nm])
    df = df[df[col_nm].str.startswith("NM_")]

    nm_to_enst: dict[str, list[str]] = defaultdict(list)
    for _, row in df.iterrows():
        nm_to_enst[row[col_nm]].append(row[col_enst])
    return dict(nm_to_enst)


def load_hcseeker_anchors(
    hcseeker_path: str,
    nm_to_enst: dict[str, list[str]],
) -> tuple[dict[str, set[int]], dict[str, set[int]]]:
    """Load HCSeeker hotspots (positive) and coldspots (negative).

    HCSeeker uses NM_ transcript IDs; we map them to ENST via nm_to_enst.

    Returns:
        (pos_anchors, neg_anchors)  each is  dict { ENST_id: set of residue positions }
    """
    logger = logging.getLogger("mcs")
    logger.info(f"Loading HCSeeker: {hcseeker_path}")

    df = pd.read_csv(hcseeker_path, sep="\t")
    logger.info(f"  {len(df)} spots ({df['type'].value_counts().to_dict()})")

    pos_anchors: dict[str, set[int]] = defaultdict(set)
    neg_anchors: dict[str, set[int]] = defaultdict(set)

    mapped, unmapped = 0, 0
    for _, row in df.iterrows():
        spot_type = row["type"]      # "hotspot" or "coldspot"
        nm_full = row["iso"]          # e.g. NM_024408.4
        nm_bare = nm_full.split(".")[0]  # strip version
        aa_start = int(row["aa_start_pos"])
        aa_end = int(row["aa_end_pos"])
        residues = set(range(aa_start, aa_end + 1))

        enst_list = nm_to_enst.get(nm_bare, [])
        if not enst_list:
            unmapped += 1
            continue
        mapped += 1

        target = pos_anchors if spot_type == "hotspot" else neg_anchors
        for enst in enst_list:
            target[enst].update(residues)

    logger.info(f"  Mapped {mapped} spots to ENST, {unmapped} unmapped (NM not in table)")
    logger.info(
        f"  Positive anchor transcripts: {len(pos_anchors)}, "
        f"Negative anchor transcripts: {len(neg_anchors)}"
    )
    return dict(pos_anchors), dict(neg_anchors)


def merge_anchors(
    rmc_pos: dict[str, set[int]],
    rmc_neg: dict[str, set[int]],
    hc_pos: dict[str, set[int]],
    hc_neg: dict[str, set[int]],
) -> tuple[dict[str, set[int]], dict[str, set[int]]]:
    """Merge RMC and HCSeeker anchors.  Remove any residues that appear in both P and N."""
    logger = logging.getLogger("mcs")

    # Union across sources
    all_tx = set(rmc_pos) | set(hc_pos) | set(rmc_neg) | set(hc_neg)
    pos_merged: dict[str, set[int]] = {}
    neg_merged: dict[str, set[int]] = {}
    conflict_count = 0

    for tx in all_tx:
        p = rmc_pos.get(tx, set()) | hc_pos.get(tx, set())
        n = rmc_neg.get(tx, set()) | hc_neg.get(tx, set())
        # Remove conflicts (residues in both positive and negative)
        overlap = p & n
        if overlap:
            conflict_count += len(overlap)
            p -= overlap
            n -= overlap
        if p:
            pos_merged[tx] = p
        if n:
            neg_merged[tx] = n

    logger.info(
        f"Merged anchors: {len(pos_merged)} pos transcripts, {len(neg_merged)} neg transcripts, "
        f"{conflict_count} conflict residues removed"
    )
    return pos_merged, neg_merged


# ═══════════════════════════════════════════════════════════════════════════════
#  Section 3: KDE smoothing (vectorised)
# ═══════════════════════════════════════════════════════════════════════════════

def kde_smooth_transcript(
    positions: np.ndarray,
    scores: np.ndarray,
    bandwidth: float,
) -> tuple[np.ndarray, np.ndarray]:
    """Gaussian KDE smoothing of composite scores along a protein.

    Args:
        positions:  sorted array of residue positions that have variants
        scores:     composite score at each position (same length as positions)
        bandwidth:  Gaussian kernel σ in residue units

    Returns:
        (query_positions, smoothed_scores)
        query_positions = every integer from min(positions) to max(positions)
    """
    if len(positions) < 2:
        return positions, scores

    query_pos = np.arange(positions[0], positions[-1] + 1)

    # sklearn KernelDensity with sample_weight gives weighted density estimate:
    #   density(x) ∝ Σ_i  weight_i · K((x - x_i)/h)
    kde = KernelDensity(kernel="gaussian", bandwidth=bandwidth)
    kde.fit(positions.reshape(-1, 1), sample_weight=scores)
    # Suppress sklearn KDE "divide by zero in log" warnings at density tails
    with warnings.catch_warnings():
        warnings.filterwarnings(
            "ignore",
            message="divide by zero",
            category=RuntimeWarning,
        )
        log_density = kde.score_samples(query_pos.reshape(-1, 1))
    smoothed = np.exp(log_density)

    return query_pos, smoothed


def smooth_all_transcripts(
    composite_data: dict[str, tuple[np.ndarray, np.ndarray, dict]],
    bandwidth: float,
    n_processes: int = 4,
) -> dict[str, tuple[np.ndarray, np.ndarray]]:
    """Smooth every transcript in parallel.

    Returns:
        { transcript_id: (query_positions, smoothed_scores) }
    """
    # Build argument list
    tx_ids = list(composite_data.keys())
    args_list = [
        (composite_data[tx][0], composite_data[tx][1], bandwidth)
        for tx in tx_ids
    ]

    with mp.Pool(n_processes) as pool:
        results = pool.starmap(kde_smooth_transcript, args_list)

    return {tx: res for tx, res in zip(tx_ids, results)}


# ═══════════════════════════════════════════════════════════════════════════════
#  Section 4: Anchor score collection & threshold utilities
# ═══════════════════════════════════════════════════════════════════════════════

def collect_anchor_scores(
    smoothed: dict[str, tuple[np.ndarray, np.ndarray]],
    pos_anchors: dict[str, set[int]],
    neg_anchors: dict[str, set[int]],
    transcripts: set[str] | None = None,
) -> tuple[np.ndarray, np.ndarray]:
    """Collect smoothed scores at positive and negative anchor residues.

    Args:
        transcripts: If None, use all transcripts in smoothed.

    Returns:
        (all_pos_scores, all_neg_scores)  as 1D arrays
    """
    pos_scores_list = []
    neg_scores_list = []

    for tx, (qpos, sscores) in smoothed.items():
        if transcripts is not None and tx not in transcripts:
            continue
        pos_min = qpos[0]
        # Positive anchors
        if tx in pos_anchors:
            anchor_positions = np.array(sorted(pos_anchors[tx]))
            valid = (anchor_positions >= qpos[0]) & (anchor_positions <= qpos[-1])
            if valid.any():
                indices = anchor_positions[valid] - pos_min
                pos_scores_list.append(sscores[indices])
        # Negative anchors
        if tx in neg_anchors:
            anchor_positions = np.array(sorted(neg_anchors[tx]))
            valid = (anchor_positions >= qpos[0]) & (anchor_positions <= qpos[-1])
            if valid.any():
                indices = anchor_positions[valid] - pos_min
                neg_scores_list.append(sscores[indices])

    all_pos = np.concatenate(pos_scores_list) if pos_scores_list else np.array([])
    all_neg = np.concatenate(neg_scores_list) if neg_scores_list else np.array([])
    return all_pos, all_neg


def dual_percentile_threshold(
    smoothed: dict[str, tuple[np.ndarray, np.ndarray]],
    pos_anchors: dict[str, set[int]],
    neg_anchors: dict[str, set[int]],
    alpha: float = 0.05,
) -> tuple[float, float, float]:
    """Compute threshold from anchor-region score distributions.

    t_pos = bottom-α percentile of smoothed scores at positive-anchor residues
    t_neg = top-α    percentile of smoothed scores at negative-anchor residues
    t_h   = max(t_pos, t_neg)

    Returns:
        (t_h, t_pos, t_neg)
    """
    all_pos, all_neg = collect_anchor_scores(smoothed, pos_anchors, neg_anchors)

    if len(all_pos) == 0:
        raise ValueError("No positive-anchor residues found in smoothed data — cannot compute threshold")
    if len(all_neg) == 0:
        raise ValueError("No negative-anchor residues found in smoothed data — cannot compute threshold")

    t_pos = np.percentile(all_pos, alpha * 100)        # bottom α of positives
    t_neg = np.percentile(all_neg, (1 - alpha) * 100)  # top α of negatives
    t_h = max(t_pos, t_neg)

    return t_h, t_pos, t_neg


# ═══════════════════════════════════════════════════════════════════════════════
#  Section 5: Objective function (Precision − η · Spillover)
# ═══════════════════════════════════════════════════════════════════════════════

def evaluate_objective(
    smoothed: dict[str, tuple[np.ndarray, np.ndarray]],
    threshold: float,
    pos_anchors: dict[str, set[int]],
    neg_anchors: dict[str, set[int]],
    eta: float = 1.0,
    transcripts: set[str] | None = None,
) -> tuple[float, float, float]:
    """Compute Q = mean_Precision - η·mean_Spillover.

    Precision is averaged only across transcripts that possess at least one
    positive anchor.  Spillover C is averaged only across transcripts that
    possess at least one negative anchor.  Transcripts lacking anchors of a
    given polarity are excluded from that metric so they cannot dilute or
    inflate the estimate.

    Precision_g = |A_g ∩ P_g| / |A_g|  (fraction of predicted that are positive anchors)
    C_g         = |A_g ∩ N_g| / |N_g|  (spillover into negative anchors)

    Args:
        transcripts: If None, evaluate over all transcripts in smoothed.

    Returns:
        (Q_mean, Prec_mean, C_mean)
    """
    if transcripts is None:
        transcripts = set(smoothed.keys())

    prec_values = []  # only transcripts with ≥1 positive anchor
    c_values = []     # only transcripts with ≥1 negative anchor

    for tx in transcripts:
        if tx not in smoothed:
            continue
        qpos, sscores = smoothed[tx]

        # Predicted constrained residues (above threshold)
        predicted = set(qpos[sscores > threshold])

        # Precision: only for transcripts that have positive anchors
        p_g = pos_anchors.get(tx, set())
        if p_g:  # transcript has ≥1 positive anchor
            if predicted:
                prec_g = len(predicted & p_g) / len(predicted)
            else:
                prec_g = 0.0  # no prediction → zero precision
            prec_values.append(prec_g)

        # Spillover: only for transcripts that have negative anchors
        n_g = neg_anchors.get(tx, set())
        if n_g:  # transcript has ≥1 negative anchor
            c_g = len(predicted & n_g) / len(n_g)
            c_values.append(c_g)

    mean_prec = float(np.mean(prec_values)) if prec_values else 0.0
    mean_c = float(np.mean(c_values)) if c_values else 0.0
    q = mean_prec - eta * mean_c

    return q, mean_prec, mean_c


# ═══════════════════════════════════════════════════════════════════════════════
#  Section 6: 2D Cross-validation — joint (bandwidth, threshold) selection
# ═══════════════════════════════════════════════════════════════════════════════

def _evaluate_single_bandwidth(
    h: float,
    composite_data: dict[str, tuple[np.ndarray, np.ndarray, dict]],
    pos_anchors: dict[str, set[int]],
    neg_anchors: dict[str, set[int]],
    folds: list[list[str]],
    cv_txs: list[str],
    eta: float,
    n_inner_proc: int,
    n_threshold_candidates: int,
) -> dict:
    """Evaluate a single bandwidth across all CV folds with threshold grid search.

    For each fold:
      1. Smooth all transcripts with bandwidth h.
      2. Collect anchor scores from TRAIN fold → determine threshold range [t_lo, t_hi].
         t_lo = 1%-percentile of positive-anchor scores (very permissive floor)
         t_hi = 95%-percentile of negative-anchor scores (stringent ceiling)
      3. Grid-search over n_threshold_candidates thresholds in [t_lo, t_hi].
      4. Pick best threshold per fold, evaluate on TEST fold.

    Per-threshold diagnostics on TEST anchors additionally include:
      - TPR_anchor(t), FPR_anchor(t), and AUC_thr(t) using binary predictions
        y_pred = 1(smoothed_score > t).
      - For a binary predictor, ROC AUC equals balanced accuracy:
            AUC_thr = (TPR + TNR) / 2
        which is directly interpretable as class-balanced discrimination.

    Returns dict with per-fold and aggregated results.
    """
    logger = logging.getLogger("mcs")
    t0 = time.time()

    # Smooth ALL transcripts once for this bandwidth
    smoothed = smooth_all_transcripts(composite_data, bandwidth=h, n_processes=n_inner_proc)

    fold_results = []
    # Aggregate best-threshold Q/J/C across folds
    q_scores, prec_scores, c_scores = [], [], []
    best_t_per_fold = []

    for fold_idx, test_txs in enumerate(folds):
        train_txs = set(cv_txs) - set(test_txs)
        test_tx_set = set(test_txs)

        # Collect anchor scores from TRAIN transcripts
        train_smoothed = {tx: smoothed[tx] for tx in train_txs if tx in smoothed}
        all_pos, all_neg = collect_anchor_scores(train_smoothed, pos_anchors, neg_anchors)

        if len(all_pos) == 0 or len(all_neg) == 0:
            logger.warning(f"  h={h} Fold {fold_idx}: insufficient anchor scores — skipping")
            continue

        # Threshold range: 1% of pos (low end) to 95% of neg (high end)
        t_lo = float(np.percentile(all_pos, 1.0))
        t_hi = float(np.percentile(all_neg, 95.0))
        # Ensure t_lo < t_hi; if not, use the narrower valid range
        if t_lo >= t_hi:
            t_lo, t_hi = min(t_lo, t_hi), max(t_lo, t_hi)
            if t_lo == t_hi:
                t_hi = t_lo * 1.1 + 1e-10  # tiny offset to avoid degenerate grid

        t_candidates = np.linspace(t_lo, t_hi, n_threshold_candidates)

        # Evaluate each threshold candidate on the TEST fold
        test_smoothed = {tx: smoothed[tx] for tx in test_tx_set if tx in smoothed}
        test_pos_scores, test_neg_scores = collect_anchor_scores(test_smoothed, pos_anchors, neg_anchors)
        y_true_test = np.concatenate(
            [
                np.ones(len(test_pos_scores), dtype=np.int8),
                np.zeros(len(test_neg_scores), dtype=np.int8),
            ]
        )
        test_scores = np.concatenate([test_pos_scores, test_neg_scores])

        all_t_results = []
        fold_best_q, fold_best_t = -np.inf, t_candidates[0]
        fold_best_prec, fold_best_c = 0.0, 0.0

        for t in t_candidates:
            q_val, prec_val, c_val = evaluate_objective(
                smoothed, t, pos_anchors, neg_anchors, eta=eta, transcripts=test_tx_set
            )
            # Compute quantile of t in both anchor distributions
            pq = float(percentileofscore(all_pos, t, kind="rank"))  # % of pos scores ≤ t
            nq = float(percentileofscore(all_neg, t, kind="rank"))  # % of neg scores ≤ t

            # TEST-anchor thresholded classification metrics
            # NOTE: ROC AUC computed on binary predictions equals balanced accuracy.
            # This turns AUC_thr into an interpretable class-balanced performance metric
            # for each fixed threshold candidate t.
            if y_true_test.size > 0:
                y_pred_bin = (test_scores > t).astype(np.int8)
                pos_mask = y_true_test == 1
                neg_mask = ~pos_mask
                tpr_anchor = float(y_pred_bin[pos_mask].mean()) if pos_mask.any() else float("nan")
                fpr_anchor = float(y_pred_bin[neg_mask].mean()) if neg_mask.any() else float("nan")
                if pos_mask.any() and neg_mask.any():
                    auc_thr = float(roc_auc_score(y_true_test, y_pred_bin))
                else:
                    auc_thr = float("nan")
            else:
                tpr_anchor, fpr_anchor, auc_thr = float("nan"), float("nan"), float("nan")

            all_t_results.append(
                (
                    float(t),
                    q_val,
                    prec_val,
                    c_val,
                    pq,
                    nq,
                    tpr_anchor,
                    fpr_anchor,
                    auc_thr,
                )
            )

            if q_val > fold_best_q:
                fold_best_q = q_val
                fold_best_t = float(t)
                fold_best_prec = prec_val
                fold_best_c = c_val

        # Quantiles of the best threshold
        best_pq = float(percentileofscore(all_pos, fold_best_t, kind="rank"))
        best_nq = float(percentileofscore(all_neg, fold_best_t, kind="rank"))

        fold_results.append({
            "fold": fold_idx,
            "t_range": (t_lo, t_hi),
            "best_t": fold_best_t,
            "best_t_pos_quantile": best_pq,
            "best_t_neg_quantile": best_nq,
            "Q": fold_best_q, "Prec": fold_best_prec, "C": fold_best_c,
            "n_test": len(test_txs),
            "all_t_results": all_t_results,
        })
        q_scores.append(fold_best_q)
        prec_scores.append(fold_best_prec)
        c_scores.append(fold_best_c)
        best_t_per_fold.append(fold_best_t)

        logger.debug(
            f"  h={h} Fold {fold_idx}: best_t={fold_best_t:.6f} "
            f"(pos_q={best_pq:.1f}%, neg_q={best_nq:.1f}%) | "
            f"Q={fold_best_q:.4f} Prec={fold_best_prec:.4f} C={fold_best_c:.4f} | "
            f"t_range=[{t_lo:.6f}, {t_hi:.6f}]"
        )

    elapsed = time.time() - t0

    if not q_scores:
        return {"h": h, "valid": False, "elapsed": elapsed}

    q_cv = float(np.mean(q_scores))
    prec_cv = float(np.mean(prec_scores))
    c_cv = float(np.mean(c_scores))
    t_median = float(np.median(best_t_per_fold))

    # Compute quantile of median threshold in full anchor distributions
    full_pos, full_neg = collect_anchor_scores(smoothed, pos_anchors, neg_anchors)
    t_med_pq = float(percentileofscore(full_pos, t_median, kind="rank")) if len(full_pos) > 0 else float("nan")
    t_med_nq = float(percentileofscore(full_neg, t_median, kind="rank")) if len(full_neg) > 0 else float("nan")
    if len(full_pos) > 0 and len(full_neg) > 0:
        y_true_full = np.concatenate([np.ones(len(full_pos), dtype=np.int8), np.zeros(len(full_neg), dtype=np.int8)])
        y_score_full = np.concatenate([full_pos, full_neg])
        anchor_auc_continuous = float(roc_auc_score(y_true_full, y_score_full))
    else:
        anchor_auc_continuous = float("nan")

    logger.info(
        f"  h={h}: best_t_median={t_median:.6f} (pos_q={t_med_pq:.1f}%, neg_q={t_med_nq:.1f}%) | "
        f"Q_CV={q_cv:.4f} (Prec={prec_cv:.4f}, C={c_cv:.4f}, anchor_auc_cont={anchor_auc_continuous:.4f}) [{elapsed:.1f}s]"
    )

    return {
        "h": h,
        "valid": True,
        "Q_CV": q_cv, "Prec_CV": prec_cv, "C_CV": c_cv,
        "anchor_auc_continuous": anchor_auc_continuous,
        "best_t_median": t_median,
        "best_t_pos_quantile": t_med_pq,
        "best_t_neg_quantile": t_med_nq,
        "per_fold": fold_results,
        "elapsed": elapsed,
    }


def cross_validate_bandwidth(
    composite_data: dict[str, tuple[np.ndarray, np.ndarray, dict]],
    pos_anchors: dict[str, set[int]],
    neg_anchors: dict[str, set[int]],
    bandwidth_candidates: list[float],
    cv_folds: int = 5,
    eta: float = 1.0,
    n_processes: int = 4,
    n_bandwidth_parallel: int = 4,
    n_threshold_candidates: int = 20,
) -> tuple[float, float, dict]:
    """Select optimal (bandwidth, threshold) via K-fold CV with 2D grid search.

    Bandwidths are evaluated in parallel using ProcessPoolExecutor.
    For each bandwidth, threshold candidates are grid-searched within each fold.

    Returns:
        (h_star, t_star, diagnostics_dict)
    """
    logger = logging.getLogger("mcs")

    # Filter to transcripts with at least one anchor
    anchor_txs = sorted(set(pos_anchors.keys()) | set(neg_anchors.keys()))
    cv_txs = [tx for tx in anchor_txs if tx in composite_data]
    logger.info(f"CV: {len(cv_txs)} transcripts with anchors (of {len(composite_data)} total)")
    logger.info(
        f"  Grid: {len(bandwidth_candidates)} bandwidths × {n_threshold_candidates} thresholds, "
        f"{n_bandwidth_parallel} parallel bandwidth workers, "
        f"{n_processes // n_bandwidth_parallel} inner processes each"
    )

    if len(cv_txs) < cv_folds:
        logger.warning(f"Fewer anchor transcripts ({len(cv_txs)}) than folds ({cv_folds}); reducing folds")
        cv_folds = max(2, len(cv_txs))

    # Assign folds (deterministic shuffle for reproducibility)
    rng = np.random.RandomState(42)
    fold_indices = rng.permutation(len(cv_txs))
    folds = [[] for _ in range(cv_folds)]
    for i, idx in enumerate(fold_indices):
        folds[i % cv_folds].append(cv_txs[idx])

    # Inner process count per bandwidth worker
    n_inner = max(1, n_processes // n_bandwidth_parallel)

    diagnostics = {
        "bandwidth_candidates": bandwidth_candidates,
        "cv_folds": cv_folds,
        "n_threshold_candidates": n_threshold_candidates,
        "results": {},
    }

    best_h, best_t, best_q = None, None, -np.inf

    # Parallel bandwidth evaluation
    logger.info(f"Launching {len(bandwidth_candidates)} bandwidth evaluations "
                f"({n_bandwidth_parallel} in parallel, {n_inner} inner proc each)")

    with ProcessPoolExecutor(max_workers=n_bandwidth_parallel) as executor:
        futures = {
            executor.submit(
                _evaluate_single_bandwidth,
                h, composite_data, pos_anchors, neg_anchors,
                folds, cv_txs, eta, n_inner, n_threshold_candidates,
            ): h
            for h in bandwidth_candidates
        }

        for future in as_completed(futures):
            h = futures[future]
            try:
                result = future.result()
            except Exception as e:
                logger.error(f"  h={h}: failed with {e}")
                continue

            if not result.get("valid", False):
                logger.warning(f"  h={h}: no valid folds — skipping")
                continue

            diagnostics["results"][h] = result

            if result["Q_CV"] > best_q:
                best_q = result["Q_CV"]
                best_h = h
                best_t = result["best_t_median"]

    # Log summary table sorted by bandwidth
    logger.info("─" * 80)
    logger.info(f"{'h':>5s} │ {'Q_CV':>8s} │ {'Prec_CV':>8s} │ {'C_CV':>8s} │ {'best_t':>10s} │ {'pos_q%':>7s} │ {'neg_q%':>7s}")
    logger.info("─" * 80)
    for h in sorted(diagnostics["results"].keys()):
        r = diagnostics["results"][h]
        marker = " ◄" if h == best_h else ""
        logger.info(
            f"{h:5.1f} │ {r['Q_CV']:8.4f} │ {r['Prec_CV']:8.4f} │ {r['C_CV']:8.4f} │ "
            f"{r['best_t_median']:10.6f} │ {r['best_t_pos_quantile']:6.1f}% │ "
            f"{r['best_t_neg_quantile']:6.1f}%{marker}"
        )
    logger.info("─" * 80)

    logger.info(f"Selected h*={best_h}, t*={best_t:.6f} with Q_CV={best_q:.4f}")
    diagnostics["h_star"] = best_h
    diagnostics["t_star"] = best_t
    diagnostics["best_Q_CV"] = best_q

    return best_h, best_t, diagnostics


# ═══════════════════════════════════════════════════════════════════════════════
#  Section 7: Final segment extraction
# ═══════════════════════════════════════════════════════════════════════════════

def extract_segments(
    composite_data: dict[str, tuple[np.ndarray, np.ndarray, dict]],
    pos_anchors: dict[str, set[int]],
    neg_anchors: dict[str, set[int]],
    bandwidth: float,
    threshold: float | None = None,
    alpha: float = 0.05,
    min_segment_len: int = 3,
    n_processes: int = 4,
) -> dict[str, set[str]]:
    """Extract constrained segments for all transcripts using given bandwidth.

    Args:
        threshold: If provided, use this threshold directly. Otherwise compute
                   via dual-percentile from anchor regions.

    Returns:
        { transcript_id: {"L91", "A92", ...} }
    """
    logger = logging.getLogger("mcs")
    logger.info(f"Extracting segments with h*={bandwidth}, min_len={min_segment_len}")

    # Smooth all transcripts
    smoothed = smooth_all_transcripts(composite_data, bandwidth=bandwidth, n_processes=n_processes)

    # Determine threshold
    if threshold is not None:
        t_h = threshold
        # Log quantile info for the explicit threshold
        all_pos, all_neg = collect_anchor_scores(smoothed, pos_anchors, neg_anchors)
        if len(all_pos) > 0 and len(all_neg) > 0:
            pq = float(percentileofscore(all_pos, t_h, kind="rank"))
            nq = float(percentileofscore(all_neg, t_h, kind="rank"))
            logger.info(f"  Using CV-selected threshold: t*={t_h:.6f} (pos_q={pq:.1f}%, neg_q={nq:.1f}%)")
        else:
            logger.info(f"  Using CV-selected threshold: t*={t_h:.6f}")
    else:
        t_h, t_pos, t_neg = dual_percentile_threshold(
            smoothed, pos_anchors, neg_anchors, alpha=alpha
        )
        logger.info(f"  Dual-percentile threshold: t_h={t_h:.6f} (t_pos={t_pos:.6f}, t_neg={t_neg:.6f})")

    # Extract contiguous segments per transcript
    results = {}
    total_segments = 0
    total_residues = 0

    for tx, (qpos, sscores) in smoothed.items():
        ref_aa_map = composite_data[tx][2]  # {pos: ref_aa}

        # Binary mask: which query positions exceed threshold
        above = sscores > t_h

        # Find contiguous runs of True values (vectorised)
        if not above.any():
            continue

        # Diff-based run detection: transitions mark segment boundaries
        padded = np.concatenate([[False], above, [False]])
        diffs = np.diff(padded.astype(int))
        starts = np.where(diffs == 1)[0]   # segment start indices
        ends = np.where(diffs == -1)[0]    # segment end indices (exclusive)

        segment_residues = set()
        seg_count = 0
        for s, e in zip(starts, ends):
            seg_positions = qpos[s:e]
            if len(seg_positions) < min_segment_len:
                continue
            seg_count += 1
            for p in seg_positions:
                if int(p) in ref_aa_map:
                    segment_residues.add(f"{ref_aa_map[int(p)]}{p}")

        if segment_residues:
            results[tx] = segment_residues
            total_segments += seg_count
            total_residues += len(segment_residues)

    logger.info(
        f"  Extracted {total_segments:,} segments across {len(results):,} transcripts "
        f"({total_residues:,} residues total)"
    )
    return results


# ═══════════════════════════════════════════════════════════════════════════════
#  Section 8: CLI
# ═══════════════════════════════════════════════════════════════════════════════

def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Missense-constrained segment extraction via anchor-guided cross-validated KDE.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    # Required
    parser.add_argument("--vcf_path", required=True, help="AlphaMissense VEP-annotated VCF (.vcf.gz)")
    parser.add_argument("--output", required=True, help="Output pickle path (.pkl)")

    # Anchor data
    parser.add_argument("--rmc_tsv", required=True, help="gnomAD RMC flat TSV (hg19 or hg38)")
    parser.add_argument("--hcseeker_tsv", required=True, help="HCSeeker spots TSV (hg19 or hg38)")
    parser.add_argument(
        "--nm_to_enst_tsv",
        default="/paedyl01/disk1/yangyxt/public_data/gene_annotation/refseq_ensembl_convert_table.txt",
        help="RefSeq NM → Ensembl ENST mapping table",
    )

    # CV parameters
    parser.add_argument(
        "--bandwidth_candidates", default="1,2,3,4,5,6,7,8,9,10,12,15,18,20,25,30",
        help="Comma-separated candidate bandwidths (AA residue units)",
    )
    parser.add_argument("--cv_folds", type=int, default=5, help="Number of CV folds")
    parser.add_argument("--eta", type=float, default=1.0, help="Negative anchor penalty weight")
    parser.add_argument("--alpha", type=float, default=0.05, help="Confidence level for dual-percentile threshold (fallback)")
    parser.add_argument("--n_threshold_candidates", type=int, default=20,
                        help="Number of threshold candidates in [1%%-pos, 99%%-neg] range")
    parser.add_argument("--n_bandwidth_parallel", type=int, default=4,
                        help="Number of bandwidths evaluated in parallel")

    # RMC classification thresholds
    parser.add_argument("--oe_constrained", type=float, default=0.4, help="RMC oe upper bound for constrained")
    parser.add_argument("--oe_unconstrained", type=float, default=0.8, help="RMC oe lower bound for unconstrained")
    parser.add_argument("--p_bh_threshold", type=float, default=0.05, help="BH-corrected p-value threshold")

    # Segment extraction
    parser.add_argument("--min_segment_len", type=int, default=3, help="Minimum contiguous segment length")

    # Runtime
    parser.add_argument("--processes", type=int, default=4, help="Number of parallel workers")
    parser.add_argument("--log_file", default=None, help="Path to log file (real-time logging)")
    parser.add_argument("--diagnostics_pkl", default=None, help="Path to save CV diagnostics pickle")

    # Diagnostic modes
    parser.add_argument("--feature_audit", action="store_true",
                        help="Compute AUC separation for candidate per-residue features, then exit (no CV/segments)")

    return parser


def main():
    parser = build_parser()
    args = parser.parse_args()

    # Setup logging (console + optional file)
    log_file = args.log_file
    if log_file is None:
        # Default: place log next to output
        log_file = str(Path(args.output).with_suffix(".log"))
    logger = setup_logging(log_file)

    logger.info("=" * 72)
    logger.info("Missense-Constrained Segment Extraction")
    logger.info("=" * 72)
    logger.info(f"Arguments: {vars(args)}")
    t_start = time.time()

    # ── Step 1: Parse VCF ──
    logger.info("─" * 40 + " Step 1: Parse VCF " + "─" * 40)
    transcript_variants = parse_vcf(args.vcf_path, n_processes=args.processes)

    # ── Step 3 (loaded early for audit path): Load anchors ──
    logger.info("─" * 40 + " Step 3: Load anchors " + "─" * 37)
    rmc_pos, rmc_neg = load_rmc_anchors(
        args.rmc_tsv,
        oe_constrained=args.oe_constrained,
        oe_unconstrained=args.oe_unconstrained,
        p_bh_threshold=args.p_bh_threshold,
    )
    nm_to_enst = _load_nm_to_enst_map(args.nm_to_enst_tsv)
    hc_pos, hc_neg = load_hcseeker_anchors(args.hcseeker_tsv, nm_to_enst)
    pos_anchors, neg_anchors = merge_anchors(rmc_pos, rmc_neg, hc_pos, hc_neg)

    # ── Feature audit mode: compute extended features, measure separation, exit ──
    if args.feature_audit:
        logger.info("─" * 40 + " Feature Audit " + "─" * 44)
        extended_data = compute_extended_features(transcript_variants)
        audit_feature_separation(extended_data, pos_anchors, neg_anchors)
        elapsed = time.time() - t_start
        logger.info(f"Feature audit complete. Total time: {elapsed:.1f}s ({elapsed / 60:.1f}min)")
        return

    # ── Step 2: Compute composite scores ──
    logger.info("─" * 40 + " Step 2: Composite scores " + "─" * 34)
    composite_data = compute_composite_scores(transcript_variants)
    del transcript_variants  # free memory

    # ── Step 4: Cross-validation bandwidth + threshold selection ──
    logger.info("─" * 40 + " Step 4: CV bandwidth + threshold " + "─" * 33)
    bandwidth_candidates = [float(x) for x in args.bandwidth_candidates.split(",")]
    h_star, t_star, diagnostics = cross_validate_bandwidth(
        composite_data, pos_anchors, neg_anchors,
        bandwidth_candidates=bandwidth_candidates,
        cv_folds=args.cv_folds,
        eta=args.eta,
        n_processes=args.processes,
        n_bandwidth_parallel=args.n_bandwidth_parallel,
        n_threshold_candidates=args.n_threshold_candidates,
    )

    if args.diagnostics_pkl:
        with open(args.diagnostics_pkl, "wb") as f:
            pickle.dump(diagnostics, f)
        logger.info(f"CV diagnostics saved to {args.diagnostics_pkl}")

    # ── Step 5: Extract final segments ──
    logger.info("─" * 40 + " Step 5: Segment extraction " + "─" * 32)
    results = extract_segments(
        composite_data, pos_anchors, neg_anchors,
        bandwidth=h_star,
        threshold=t_star,
        alpha=args.alpha,
        min_segment_len=args.min_segment_len,
        n_processes=args.processes,
    )

    # ── Save output ──
    with open(args.output, "wb") as f:
        pickle.dump(results, f)

    elapsed = time.time() - t_start
    logger.info("=" * 72)
    logger.info(f"Done. {len(results):,} transcripts with segments → {args.output}")
    logger.info(f"Total time: {elapsed:.1f}s ({elapsed / 60:.1f}min)")
    logger.info("=" * 72)


if __name__ == "__main__":
    main()
