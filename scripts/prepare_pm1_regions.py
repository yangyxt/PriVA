#!/usr/bin/env python3
"""Prepare all PM1 evidence regions for PriVA ACMG variant classification.

PM1 has three arms:
  1. DAS-selected intolerant domains — InterPro domains whose AlphaMissense
     score distribution is statistically indistinguishable from known
     functional domains, after excluding coldspot-overlapping, benign-ClinVar,
     and common-gnomAD domains.
  2. HCSeeker hotspots — ClinVar-based observed mutational hotspots (amino
     acid stretches enriched for pathogenic variants relative to benign).
  3. gnomAD RMC constrained regions — sub-genic regions depleted of
     missense variation in population data (observed/expected < threshold
     with FDR-corrected significance).

Output: a single pickle file containing all three region sets plus metadata.
"""

import argparse
import logging
import os
import pickle
import re
import sys
from collections import defaultdict

import yaml

logger = logging.getLogger("pm1_regions")

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
PRIVA_DIR = os.path.dirname(SCRIPT_DIR)

DEFAULT_NM_TO_ENST_TSV = (
    "/paedyl01/disk1/yangyxt/public_data/gene_annotation/"
    "refseq_ensembl_convert_table.txt"
)


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _load_nm_to_enst_map(mapping_path: str) -> dict:
    """Load RefSeq NM → Ensembl ENST mapping table.

    Parameters
    ----------
    mapping_path : str
        Path to the TSV with columns:
        Gene stable ID, Transcript stable ID, HGNC symbol,
        UniProtKB/Swiss-Prot ID, RefSeq mRNA ID.

    Returns
    -------
    dict[str, list[str]]
        {NM_id_without_version: [ENST_id, ...]}
    """
    nm_to_enst: dict = defaultdict(list)
    with open(mapping_path) as fh:
        header = fh.readline()  # skip header
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 5:
                continue
            enst = parts[1].strip()
            refseq = parts[4].strip()
            if not refseq.startswith("NM_"):
                continue
            nm_no_ver = refseq.split(".")[0]
            if enst not in nm_to_enst[nm_no_ver]:
                nm_to_enst[nm_no_ver].append(enst)
    logger.info(
        "Loaded NM→ENST mapping: %d NM accessions → %d unique ENST IDs",
        len(nm_to_enst),
        len({e for lst in nm_to_enst.values() for e in lst}),
    )
    return dict(nm_to_enst)


def _parse_rmc_aa_pos(raw: str) -> int:
    """Parse RMC amino-acid position strings like 'Arg97' → 97."""
    m = re.search(r"(\d+)", raw)
    if m is None:
        raise ValueError(f"Cannot parse AA position from '{raw}'")
    return int(m.group(1))


# ---------------------------------------------------------------------------
# HCSeeker hotspots
# ---------------------------------------------------------------------------

def load_hcseeker_hotspots(
    hcseeker_tsv: str,
    nm_to_enst: dict,
) -> dict:
    """Load HCSeeker hotspots → {ENST: set(AA positions)}.

    Only rows with type == 'hotspot' are kept.
    """
    result: dict = defaultdict(set)
    mapped = 0
    unmapped = 0

    with open(hcseeker_tsv) as fh:
        header = fh.readline()
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 7:
                continue
            spot_type = parts[0].strip()
            if spot_type != "hotspot":
                continue
            iso = parts[4].strip()
            aa_start = int(parts[5])
            aa_end = int(parts[6])

            nm_no_ver = iso.split(".")[0]
            enst_list = nm_to_enst.get(nm_no_ver)
            if not enst_list:
                unmapped += 1
                continue

            mapped += 1
            positions = set(range(aa_start, aa_end + 1))
            for enst in enst_list:
                result[enst].update(positions)

    total_residues = sum(len(s) for s in result.values())
    logger.info(
        "HCSeeker hotspots: %d mapped, %d unmapped, "
        "%d transcripts, %d residues",
        mapped, unmapped, len(result), total_residues,
    )
    return dict(result)


# ---------------------------------------------------------------------------
# gnomAD RMC constrained regions
# ---------------------------------------------------------------------------

def load_rmc_constrained(
    rmc_tsv: str,
    oe_threshold: float = 0.4,
    p_bh_threshold: float = 0.05,
) -> dict:
    """Load RMC constrained regions → {ENST: set(AA positions)}.

    Applies Benjamini-Hochberg FDR correction on raw p-values, then
    selects regions with oe < *oe_threshold* AND p_bh < *p_bh_threshold*.
    """
    from scipy.stats import false_discovery_control

    records = []
    with open(rmc_tsv) as fh:
        header = fh.readline().rstrip("\n").split("\t")
        # Determine column indices — hg19 has 13 cols, hg38 has 10
        col = {name: idx for idx, name in enumerate(header)}
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            transcript = parts[col["transcript"]].strip()
            start_raw = parts[col["start_aa_raw"]].strip()
            end_raw = parts[col["end_aa_raw"]].strip()
            oe = float(parts[col["oe"]])
            p_val = float(parts[col["p"]])

            start_pos = _parse_rmc_aa_pos(start_raw)
            end_pos = _parse_rmc_aa_pos(end_raw)

            records.append((transcript, start_pos, end_pos, oe, p_val))

    logger.info("RMC: loaded %d regions total", len(records))

    if not records:
        return {}

    # BH FDR correction
    import numpy as np
    raw_p = np.array([r[4] for r in records], dtype=np.float64)
    p_bh = false_discovery_control(raw_p, method="bh")

    result: dict = defaultdict(set)
    constrained_count = 0
    for (transcript, start_pos, end_pos, oe, _), adj_p in zip(records, p_bh):
        if oe < oe_threshold and adj_p < p_bh_threshold:
            constrained_count += 1
            result[transcript].update(range(start_pos, end_pos + 1))

    total_residues = sum(len(s) for s in result.values())
    logger.info(
        "RMC constrained: %d / %d regions (oe<%.2f, p_bh<%.3f), "
        "%d transcripts, %d residues",
        constrained_count, len(records),
        oe_threshold, p_bh_threshold,
        len(result), total_residues,
    )
    return dict(result)


# ---------------------------------------------------------------------------
# DAS intolerant-domain pipeline
# ---------------------------------------------------------------------------

def run_das_pipeline(
    config: dict,
    output_dir: str,
    assembly: str,
    threads: int,
) -> set:
    """Run (or reload) the DAS intolerant-domain pipeline.

    Returns a set of intolerant domain-path strings, e.g.
    ``"ENSG00000198835:Pfam:PF00029"``.
    """
    pkl_name = f"all_intolerant_domains.{assembly}.pkl"
    pkl_path = os.path.join(output_dir, pkl_name)

    pickle_file = config.get("alphamissense_pd_stat", "")
    am_vcf = config.get("alphamissense_vep_vcf", "")
    hcseeker_tsv = os.path.join(
        PRIVA_DIR, "data", "HCSeeker", f"HC_spots.{assembly}.tsv"
    )
    clinvar_vcf = config.get("clinvar_vep_vcf", "") or config.get("clinvar_vcf", "")

    # Check if we can skip regeneration
    if os.path.isfile(pkl_path) and os.path.isfile(pickle_file):
        if os.path.getmtime(pkl_path) > os.path.getmtime(pickle_file):
            logger.info("DAS: loading existing %s (newer than input)", pkl_path)
            with open(pkl_path, "rb") as fh:
                return pickle.load(fh)

    # Import the DAS module from the same directory
    if SCRIPT_DIR not in sys.path:
        sys.path.insert(0, SCRIPT_DIR)
    from am_pick_intolerant_domains import analyze_domain_data

    logger.info("DAS: running intolerant-domain pipeline …")
    analyze_domain_data(
        pickle_file,
        output_dir,
        threads,
        assembly=assembly,
        am_vcf=am_vcf,
        hcseeker_spots_tsv=hcseeker_tsv,
        clinvar_vcf=clinvar_vcf,
    )

    with open(pkl_path, "rb") as fh:
        domains = pickle.load(fh)
    logger.info("DAS: %d intolerant domains", len(domains))
    return domains


# ---------------------------------------------------------------------------
# Main entry-point
# ---------------------------------------------------------------------------

def prepare_pm1_regions(args):
    """Build the combined PM1 regions pickle."""
    with open(args.config) as fh:
        config = yaml.safe_load(fh)

    assembly = args.assembly
    output_path = args.output

    # --- Paths ---
    hcseeker_tsv = os.path.join(
        PRIVA_DIR, "data", "HCSeeker", f"HC_spots.{assembly}.tsv"
    )
    rmc_tsv = os.path.join(
        PRIVA_DIR, "data", "gnomAD_RMC",
        f"gnomad_v2.1.1_rmc_flat.{assembly}.tsv",
    )
    nm_to_enst_tsv = args.nm_to_enst_tsv

    # --- NM→ENST mapping ---
    nm_to_enst = _load_nm_to_enst_map(nm_to_enst_tsv)

    # --- Arm 1: DAS intolerant domains ---
    if args.skip_das:
        das_pkl = config.get(
            "all_intolerant_domains",
            os.path.join(
                os.path.dirname(config.get("alphamissense_pd_stat", "")),
                f"all_intolerant_domains.{assembly}.pkl",
            ),
        )
        if not os.path.isfile(das_pkl):
            logger.error("--skip_das set but %s not found", das_pkl)
            sys.exit(1)
        logger.info("DAS: loading pre-computed %s", das_pkl)
        with open(das_pkl, "rb") as fh:
            das_intolerant_domains = pickle.load(fh)
    else:
        output_dir = os.path.dirname(
            config.get("alphamissense_pd_stat", output_path)
        )
        das_intolerant_domains = run_das_pipeline(
            config, output_dir, assembly, args.threads
        )

    # --- Arm 2: HCSeeker hotspots ---
    hcseeker_hotspots = load_hcseeker_hotspots(hcseeker_tsv, nm_to_enst)

    # --- Arm 3: gnomAD RMC constrained ---
    rmc_constrained = load_rmc_constrained(
        rmc_tsv,
        oe_threshold=args.rmc_oe_threshold,
        p_bh_threshold=args.rmc_p_bh_threshold,
    )

    # --- Bundle ---
    n_hcseeker_residues = sum(len(s) for s in hcseeker_hotspots.values())
    n_rmc_residues = sum(len(s) for s in rmc_constrained.values())

    pm1_regions = {
        "das_intolerant_domains": das_intolerant_domains,
        "hcseeker_hotspots": hcseeker_hotspots,
        "rmc_constrained": rmc_constrained,
        "metadata": {
            "assembly": assembly,
            "rmc_oe_threshold": args.rmc_oe_threshold,
            "rmc_p_bh_threshold": args.rmc_p_bh_threshold,
            "n_das_intolerant_domains": len(das_intolerant_domains),
            "n_hcseeker_hotspot_transcripts": len(hcseeker_hotspots),
            "n_hcseeker_hotspot_residues": n_hcseeker_residues,
            "n_rmc_constrained_transcripts": len(rmc_constrained),
            "n_rmc_constrained_residues": n_rmc_residues,
        },
    }

    with open(output_path, "wb") as fh:
        pickle.dump(pm1_regions, fh, protocol=pickle.HIGHEST_PROTOCOL)
    logger.info("Saved PM1 regions → %s", output_path)

    # Summary table
    meta = pm1_regions["metadata"]
    print("\n===== PM1 Region Summary =====")
    print(f"  Assembly:                        {meta['assembly']}")
    print(f"  DAS intolerant domains:          {meta['n_das_intolerant_domains']}")
    print(f"  HCSeeker hotspot transcripts:    {meta['n_hcseeker_hotspot_transcripts']}")
    print(f"  HCSeeker hotspot residues:       {meta['n_hcseeker_hotspot_residues']}")
    print(f"  RMC constrained transcripts:     {meta['n_rmc_constrained_transcripts']}")
    print(f"  RMC constrained residues:        {meta['n_rmc_constrained_residues']}")
    print(f"  RMC OE threshold:                {meta['rmc_oe_threshold']}")
    print(f"  RMC p_BH threshold:              {meta['rmc_p_bh_threshold']}")
    print("==============================\n")


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Prepare PM1 evidence regions for PriVA ACMG classification."
    )
    parser.add_argument(
        "--config", required=True,
        help="Path to PriVA config YAML.",
    )
    parser.add_argument(
        "--output", required=True,
        help="Output pickle path, e.g. pm1_regions.hg19.pkl.",
    )
    parser.add_argument(
        "--assembly", required=True, choices=["hg19", "hg38"],
        help="Genome assembly (hg19 or hg38).",
    )
    parser.add_argument(
        "--threads", type=int, default=12,
        help="Threads for DAS pipeline (default: 12).",
    )
    parser.add_argument(
        "--rmc_oe_threshold", type=float, default=0.4,
        help="RMC observed/expected threshold (default: 0.4).",
    )
    parser.add_argument(
        "--rmc_p_bh_threshold", type=float, default=0.05,
        help="RMC BH-corrected p-value threshold (default: 0.05).",
    )
    parser.add_argument(
        "--skip_das", action="store_true",
        help="Skip DAS pipeline; load existing all_intolerant_domains pkl.",
    )
    parser.add_argument(
        "--nm_to_enst_tsv", default=DEFAULT_NM_TO_ENST_TSV,
        help="Path to RefSeq NM → Ensembl ENST mapping TSV.",
    )
    return parser


if __name__ == "__main__":
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s [%(name)s] %(levelname)s: %(message)s",
    )
    parser = build_parser()
    args = parser.parse_args()
    prepare_pm1_regions(args)
