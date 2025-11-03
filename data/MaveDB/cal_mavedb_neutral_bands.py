#!/usr/bin/env python3
"""
compute_neutral_bands.py

Scan a directory of MaveDB-style score CSV files, extract synonymous variants
(hgvs_pro ending with '='), and compute the 5% / 95% percentile neutral band
for each score set. If there are insufficient synonymous variants, fall back to
a robust non-parametric band using median ± k * MAD (k ≈ 2.439 for 90% central).

Usage:
    python compute_neutral_bands.py --scores-dir ./mavedb_data --output neutral_bands.tsv
    python compute_neutral_bands.py --scores-dir ./mavedb_data --output neutral_bands.tsv --nonmonotonic-list nonmono.txt

Output columns:
    urn         : urn:mavedb:... (experiment/score-set URN, without trailing #)
    q05         : 5th percentile (or NA)
    q95         : 95th percentile (or NA)
    n_syn_used  : number of synonyms used after QC (0 when MAD fallback is used)
    band_source : how the band was computed (synonym_percentiles_qc | mad90_all_variants | percentiles_all_variants | reason for NA)
    file_path   : path to score CSV processed
"""

import argparse
import os
import sys
import glob
import math
import logging
import numpy as np
import pandas as pd


logger = logging.getLogger()
console_handler = logging.StreamHandler()
console_handler.setLevel(logging.INFO)
formatter = logging.Formatter("%(levelname)s:%(asctime)s:%(message)s")
console_handler.setFormatter(formatter)
logger.addHandler(console_handler)
logger.setLevel(os.environ.get("LOGLEVEL", "INFO"))


# ----------------------- Helpers for scoring columns & URNs -----------------------

def choose_score_column(df):
    # Preference order
    for c in ['score']:
        if c in df.columns:
            return c
    # Fallback: first numeric column that doesn't look like an error/df column
    numeric_cols = df.select_dtypes(include=[np.number]).columns.tolist()
    for c in numeric_cols:
        cl = c.lower()
        if cl not in ('sd', 'se', 'df', 'sigma', 'sem', 'stderr', 'std'):
            return c
    return None


def get_accession_urn(df, filename):
    # Prefer the 'accession' column and strip anything after '#'
    if 'accession' in df.columns:
        val = df['accession'].dropna().astype(str)
        if len(val) > 0:
            for v in val.values:
                if v.startswith('urn:'):
                    return v.split('#')[0]
    # Fallback: derive from filename (e.g., urn-mavedb-00000001-d-1.scores.*)
    basename = os.path.basename(filename)
    name = os.path.splitext(basename)[0]
    if name.startswith('urn-'):
        parts = name.split('-', 2)
        if len(parts) >= 3:
            return parts[0] + ':' + parts[1] + ':' + parts[2]
    return name



# ----------------------- MAD-based neutral band (fallback) -----------------------

def neutral_band_mad(scores, alpha=0.10, min_n=50):
    """
    Robust neutral band from all variant scores using median ± k * MAD.
    alpha = 0.10 → central 90% band, k ≈ 2.439 (1.4826 * z_0.95).
    Returns (lo, hi, meta) or (nan, nan, meta) if insufficient data.
    """
    x = np.asarray(scores, dtype=float)
    x = x[np.isfinite(x)]
    if x.size < min_n:
        return (math.nan, math.nan, {"reason": "insufficient_n", "n": int(x.size)})

    m = np.median(x)
    mad = np.median(np.abs(x - m))
    if not np.isfinite(mad):
        return (math.nan, math.nan, {"reason": "mad_nonfinite", "n": int(x.size)})

    if mad == 0.0:
        eps = 1e-6
        return (float(m - eps), float(m + eps),
                {"reason": "zero_mad", "median": float(m), "mad": float(mad), "n": int(x.size)})

    Z = {0.20: 1.28155, 0.10: 1.64485, 0.05: 1.95996}
    z = Z.get(alpha, 1.64485)  # default 90%
    k = 1.4826 * z

    lo, hi = float(m - k * mad), float(m + k * mad)
    return (lo, hi, {"median": float(m), "mad": float(mad), "alpha": alpha, "k": float(k), "n": int(x.size)})


# -------- percentile-based fallback when MAD cannot be applied --------

def neutral_band_percentiles(scores, lo_q=5, hi_q=95, min_n=20):
    """
    Direct percentile band from all variant scores.
    Returns (lo, hi, meta) or (nan, nan, meta) if insufficient data.
    """
    x = np.asarray(scores, dtype=float)
    x = x[np.isfinite(x)]
    if x.size < min_n:
        return (math.nan, math.nan, {"reason": "insufficient_n", "n": int(x.size)})
    lo = float(np.nanpercentile(x, lo_q))
    hi = float(np.nanpercentile(x, hi_q))
    return (lo, hi, {"n": int(x.size), "lo_q": lo_q, "hi_q": hi_q})


# ----------------------- Main per-file processing -----------------------

def process_file(path, nonmono_set):
    # Load CSV (comma then tab)
    try:
        df = pd.read_csv(path, low_memory=False)
    except Exception:
        try:
            df = pd.read_csv(path, sep='\t', low_memory=False)
        except Exception as e:
            logger.warning(f"Could not read {path}: {e}")
            return None

    urn = get_accession_urn(df, path)

    if nonmono_set and urn in nonmono_set:
        return {
            'urn': urn, 'q05': 'NA', 'q95': 'NA',
            'n_syn_used': 0, 'band_source': 'nonmonotonic_list', 'file_path': path
        }

    score_col = choose_score_column(df)
    if score_col is None:
        logger.warning(f"No score column found in {path}, skipping.")
        return {
            'urn': urn, 'q05': 'NA', 'q95': 'NA',
            'n_syn_used': 0, 'band_source': 'no_score_column', 'file_path': path
        }

    df['_score'] = pd.to_numeric(df[score_col], errors='coerce')

    # -------- Synonym route (preferred) --------
    if 'hgvs_pro' in df.columns:
        syn = df[df['hgvs_pro'].astype(str).str.endswith('=') & df['_score'].notna()].copy()

        if len(syn) >= 20:
            logger.info(f"Using {len(syn)} synonymous variants for neutral band in {path}")
            q05 = float(np.nanpercentile(syn['_score'].values, 5))
            q95 = float(np.nanpercentile(syn['_score'].values, 95))
            logger.info(f"Computed neutral band percentiles: q05={q05}, q95={q95} for {path}")
            if q05 == q95 and "exp.score" in df.columns:
                # Edge case: all synonyms have the same score; try exp.score if available
                exp_scores = pd.to_numeric(syn['exp.score'], errors='coerce')
                exp_scores = exp_scores[np.isfinite(exp_scores)]
                if exp_scores.size >= 20:
                    q05 = float(np.nanpercentile(exp_scores.values, 5))
                    q95 = float(np.nanpercentile(exp_scores.values, 95))
                    logger.info(f"Recomputed neutral band from exp.score: q05={q05}, q95={q95} for {path}")

            if q05 != q95:
                return {
                    'urn': urn,
                    'q05': q05,
                    'q95': q95,
                    'n_syn_used': int(len(syn)),
                    'band_source': 'synonym_percentiles_qc',
                    'file_path': path
                }

    # -------- MAD fallback on ALL variants (noncoding/no-syn cases) --------
    logger.warning(f"Falling back to MAD-based neutral band for {path}, the table only contains {len(df)} total variants and {len(syn) if 'syn' in locals() else 0} synonyms.")
    all_scores = pd.to_numeric(df['_score'], errors='coerce')
    all_scores = all_scores[np.isfinite(all_scores)]
    if all_scores.size == 0:
        return {
            'urn': urn, 'q05': 'NA', 'q95': 'NA',
            'n_syn_used': 0, 'band_source': 'no_usable_scores', 'file_path': path
        }

    lo, hi, meta = neutral_band_mad(all_scores.values, alpha=0.10, min_n=50)
    if np.isfinite(lo) and np.isfinite(hi):
        logger.info(f"Computed MAD-based neutral band: q05={lo}, q95={hi} for {path}")
        return {
            'urn': urn,
            'q05': float(lo),
            'q95': float(hi),
            'n_syn_used': 0,  # MAD fallback uses all variants
            'band_source': 'mad90_all_variants',
            'file_path': path
        }

    # -------- NEW: Percentile fallback if MAD is not applicable --------
    logger.warning(f"MAD fallback not applicable for {path} (reason: {meta.get('reason', 'unknown')}); using raw 5/95 percentiles over all variants.")
    plo, phi, pmeta = neutral_band_percentiles(all_scores.values, lo_q=5, hi_q=95, min_n=20)
    if np.isfinite(plo) and np.isfinite(phi):
        logger.info(f"Computed percentile-based neutral band: q05={plo}, q95={phi} for {path}")
        return {
            'urn': urn,
            'q05': float(plo),
            'q95': float(phi),
            'n_syn_used': 0,
            'band_source': 'percentiles_all_variants',
            'file_path': path
        }

    # If even percentiles can't be computed, fall back to NA with reason
    reason = meta.get("reason", "mad_failed")
    if 'reason' in pmeta:
        reason = f"{reason}|percentiles_{pmeta['reason']}"
    return {
        'urn': urn, 'q05': 'NA', 'q95': 'NA',
        'n_syn_used': 0, 'band_source': f'mad_and_percentiles_na_{reason}', 'file_path': path
    }


def find_score_files(root_dir):
    # Zenodo archive typically has a 'csv' subfolder; include recursive glob just in case.
    patterns = ['**/*scores.csv', '*scores.csv']
    files = []
    for p in patterns:
        files.extend(glob.glob(os.path.join(root_dir, p), recursive=True))
    return sorted(set(files))


# ----------------------- Entry point -----------------------

def main():
    parser = argparse.ArgumentParser(description="Compute neutral bands from MaveDB score CSVs (synonym 5/95 or MAD fallback).")
    parser.add_argument('--scores-dir', required=True, help="Directory containing score CSVs (unzipped MaveDB archive).")
    parser.add_argument('--output', required=True, help="Output TSV file for neutral bands.")
    parser.add_argument('--nonmonotonic-list', required=False, help="Optional file with URNs (one per line) to mark as nonmonotonic (set NA).")
    args = parser.parse_args()

    nonmono_set = set()
    if args.nonmonotonic_list:
        with open(args.nonmonotonic_list) as fh:
            for line in fh:
                u = line.strip()
                if u:
                    nonmono_set.add(u)

    files = find_score_files(args.scores_dir)
    logger.info(f"Found {len(files)} candidate score files")
    results = []
    for f in files:
        try:
            res = process_file(f, nonmono_set)
            if res:
                results.append(res)
        except Exception as e:
            logger.error(f"Failed processing {f}: {e}")

    out_df = pd.DataFrame(results, columns=['urn','q05','q95','n_syn_used','band_source','file_path'])
    out_df.to_csv(args.output, sep='\t', index=False)
    logger.info(f"Wrote {len(out_df)} rows to {args.output}")


if __name__ == '__main__':
    main()
