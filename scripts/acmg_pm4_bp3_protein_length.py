#!/usr/bin/env python3
"""PM4 and BP3 -- the variant changes protein length without abolishing it.

Mutually exclusive readings of the same class of variant: in-frame
insertions and deletions, and stop-loss variants that extend the protein.

    PM4  the length change is not in a repetitive region, or it falls in a
         functional or intolerant domain -- so the change plausibly matters
    BP3  the length change is in a repetitive region of no known function --
         so it plausibly does not

find_overlaps_bedtools_efficient does the repeat-region intersection.
"""

import logging
import pandas as pd
import numpy as np
from typing import Tuple


logger = logging.getLogger(__name__)


def find_overlaps_bedtools_efficient(variants_df, regions_file, method = "any"):
    """
    Find variants that overlap with regions using pybedtools with optimal efficiency.
    Handles conversion between 1-based variant coordinates and 0-based BED coordinates.

    Args:
        variants_df: DataFrame with variants (using 1-based coordinates)
        regions_file: Path to BED file with regions (using 0-based coordinates)

    Returns:
        set: Set of variant IDs that overlap with regions
    """
    import pybedtools

    # Convert variants to BED format (1-based to 0-based)
    variants_bed = pybedtools.BedTool.from_dataframe(
        variants_df[['chrom', 'pos', 'variant_id']].assign(
            start=lambda x: x['pos'] - 1,  # Convert 1-based to 0-based
            end=lambda x: x['pos'] - 1 + variants_df.get('ref', 'A').str.len()  # End position
        )[['chrom', 'start', 'end', 'variant_id']]
    )

    # Load regions (already in 0-based BED format)
    regions_bed = pybedtools.BedTool(regions_file)

    # Find overlaps - this keeps only the features from variants_bed that overlap
    if method == "any":
        # Find intersections of any size
        intersect_result = variants_bed.intersect(regions_bed, wa=True, u=True)
    elif method == "all":
        # Find intersections where the variant interval is fully contained in the region interval
        intersect_result = variants_bed.intersect(regions_bed, wa=True, f=1, u=True)

    # Extract the variant IDs that had overlaps
    overlapping_variants = set()
    for feature in intersect_result:
        # The variant_id is the 4th field (index 3)
        overlapping_variants.add(feature[3])

    return overlapping_variants


def PM4_BP3_criteria(df: pd.DataFrame,
                     pvs1_criteria: np.ndarray,
                     repeat_regions_file: str,
                     loc_intol_domain: np.ndarray) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    Combined PM4 and BP3 criteria assignment for in-frame length-changing variants.

    PM4 and BP3 are mutually exclusive:
    - PM4: In-frame indels NOT in repetitive regions OR in functional/intolerant domains
           Also applies to stop-loss variants (protein extension)
    - BP3: In-frame indels IN repetitive regions WITHOUT functional significance

    Per ClinGen SVI guidelines (PMC6185798):
    - PM4 applies to in-frame deletions/insertions AND stop-loss variants (protein extension)
    - BP3 applies to in-frame deletions/insertions in repetitive regions without known function
    - Stop-loss variants are NOT null variants (not eligible for PVS1), they get PM4
    - PVS1 only applies to: nonsense, frameshift, canonical splice, start_lost (null variants)

    Workflow:
    1. First identify if variants are in repetitive regions
    2. Identify all in-frame length-changing variants
    3. Check functional domain involvement
    4. Assign PM4 or BP3 based on:
       - PM4: NOT in repeat OR in functional domain (or stop-loss)
       - BP3: IN repeat AND NOT in functional domain

    Args:
        df: DataFrame with variant annotations
        pvs1_criteria: Array of PVS1 strength values
        repeat_regions_file: Path to BED file with repetitive regions
        loc_intol_domain: Boolean array indicating variants in intolerant domains

    Returns:
        Tuple of:
        - pm4_array: PM4 criteria strength values (int array)
        - bp3_array: BP3 criteria strength values (int array)
        - in_repeat: Boolean array indicating variants in repetitive regions
    """
    # =========================================================================
    # Step 1: Identify variants in repetitive regions
    # =========================================================================
    in_repeat_vars = find_overlaps_bedtools_efficient(df, repeat_regions_file, method="all")
    in_repeat = df['variant_id'].isin(in_repeat_vars)

    # Also check for Low_complexity domains in DOMAINS annotation
    low_complexity_domain = df['DOMAINS'].fillna(".").str.contains('Low_complexity')
    in_repetitive_region = in_repeat | low_complexity_domain

    logger.info(f"PM4_BP3_criteria: {in_repetitive_region.sum()} variants in repetitive regions")

    # =========================================================================
    # Step 2: Identify all in-frame length-changing variants
    # =========================================================================
    # Standard in-frame indels from VEP consequences
    inframe_indels = (df['Consequence'].fillna("").str.contains('inframe_deletion')) | \
                     (df['Consequence'].fillna("").str.contains('inframe_insertion'))

    # Splicing-induced in-frame length changes (not frameshift)
    splicing_inframe = df["splicing_len_changing"].fillna(False) & \
                       ~df["splicing_frameshift"].fillna(False)

    # 5'UTR-induced in-frame length changes (not frameshift)
    utr_inframe = df["5UTR_len_changing"].fillna(False) & \
                  ~df["5UTR_frameshift"].fillna(False)

    # Also include vep_consq_len_changing but exclude frameshift and stop_gained
    # BP3 applies to IN-FRAME length changes only (not LoF variants)
    vep_inframe = df['vep_consq_len_changing'].fillna(False) & \
                  ~df["Consequence"].fillna("").str.contains("frameshift") & \
                  ~df["Consequence"].fillna("").str.contains("stop_gained")

    # Combined in-frame length-changing variants
    all_inframe = inframe_indels | splicing_inframe | utr_inframe | vep_inframe

    logger.info(f"PM4_BP3_criteria: {all_inframe.sum()} total in-frame length-changing variants")

    # =========================================================================
    # Step 3: Identify stop-loss variants (special case for PM4 only)
    # =========================================================================
    # Stop-loss variants get PM4 regardless of repetitive region
    # Stop-loss causes protein extension, not LoF
    stop_loss_variants = df["Consequence"].str.contains("stop_lost").fillna(False)
    logger.info(f"PM4_BP3_criteria: {stop_loss_variants.sum()} stop-loss variants (always get PM4)")

    # NMD-escaping truncating variants (frameshift / stop-gained in the last exon)
    # still translate a truncated protein, so ClinGen treats them as a length change
    # (PM4) rather than a null variant (PVS1). The `pvs1_criteria < 2` gate below
    # keeps them from double-counting when PVS1 is still present (LoF gene).
    nmd_escaping_trunc = (df['Consequence'].fillna("").str.contains("frameshift") |
                          df['Consequence'].fillna("").str.contains("stop_gained")) & \
                         (df['NMD'].fillna(".").str.contains("escaping") |
                          df['LoF_filter'].fillna(".").str.contains("END_TRUNC"))
    logger.info(f"PM4_BP3_criteria: {nmd_escaping_trunc.sum()} NMD-escaping truncating variants")

    # =========================================================================
    # Step 4: Identify variants with functional domain involvement
    # =========================================================================
    # These get PM4 even if in repetitive region
    has_functional_domain = loc_intol_domain | \
                            df["splicing_span_intol_domain"].fillna(False) | \
                            df["5UTR_span_intol_domain"].fillna(False) | \
                            df['span_functional_domains'].fillna(False)

    logger.info(f"PM4_BP3_criteria: {has_functional_domain.sum()} variants span functional domains")

    # =========================================================================
    # Step 5: Exclude variants that already have strong PVS1 (frameshift/NMD)
    # =========================================================================
    frameshift_with_pvs1 = df["Consequence"].fillna("").str.contains("frameshift") & \
                           (pvs1_criteria > 0) & ~stop_loss_variants
    strong_nmd = (pvs1_criteria >= 3) & ~stop_loss_variants
    exclude_from_pm4 = frameshift_with_pvs1 | strong_nmd

    # =========================================================================
    # Step 6: Assign PM4 and BP3 (mutually exclusive)
    # =========================================================================
    # PM4 logic:
    # - In-frame indels that are: (NOT in repeat) OR (in functional domain)
    # - AND not excluded by PVS1
    # - OR stop-loss variants (always get PM4)
    pm4_inframe_eligible = all_inframe & ~exclude_from_pm4 & \
                           (~in_repetitive_region | has_functional_domain)
    pm4_eligible = pm4_inframe_eligible | stop_loss_variants | nmd_escaping_trunc

    # BP3 logic:
    # - In-frame indels that are: IN repeat AND NOT in functional domain
    bp3_eligible = all_inframe & in_repetitive_region & ~has_functional_domain

    # Ensure mutual exclusivity: PM4 takes precedence
    # (If variant is in functional domain within repeat region, it gets PM4, not BP3)
    bp3_eligible = bp3_eligible & ~pm4_eligible

    # =========================================================================
    # Step 7: Create output arrays with strength levels
    # =========================================================================
    pm4_array = np.zeros(len(df), dtype=int)
    bp3_array = np.zeros(len(df), dtype=int)

    # PM4 is assigned at Moderate level (2) when PVS1 < 2
    pm4_array[pm4_eligible & (pvs1_criteria < 2)] = 2

    # BP3 is assigned at Supporting level (1)
    # Also exclude if PM4 was assigned (double-check mutual exclusivity)
    bp3_array[bp3_eligible & (pm4_array == 0)] = 1

    # =========================================================================
    # Step 8: Log summary statistics
    # =========================================================================
    logger.info(f"PM4_BP3_criteria: {(pm4_array > 0).sum()} variants assigned PM4")
    logger.info(f"PM4_BP3_criteria: {(bp3_array > 0).sum()} variants assigned BP3")

    # Verify mutual exclusivity
    both_assigned = (pm4_array > 0) & (bp3_array > 0)
    if both_assigned.any():
        logger.error(f"PM4_BP3_criteria: ERROR - {both_assigned.sum()} variants have both PM4 and BP3!")
    else:
        logger.info("PM4_BP3_criteria: Verified PM4 and BP3 are mutually exclusive")

    return pm4_array, bp3_array, in_repeat
