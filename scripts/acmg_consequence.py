#!/usr/bin/env python3
"""Shared foundation: read VEP's annotation and say what the variant does.

No ACMG criterion lives here. These are the per-variant facts that several
criteria need but none of them owns:

    vep_consq_interpret       whether the consequence is loss-of-function,
                              and whether it changes protein length
    truncate_fraction         how much of the protein a truncation removes
    get_variant_type          missense / nonsense / frameshift, from HGVSp
    coerce_pext_value         collapse VEP's ampersand-joined pext values to
                              one number

Read by PVS1, PS1/PM5 and the ranking step.
"""

import logging
import re
import multiprocessing as mp
import pandas as pd
import numpy as np
from typing import Tuple, Dict


logger = logging.getLogger(__name__)


def coerce_pext_value(value):
    """
    Convert VEP BigWig PEXT annotations to a single numeric value.

    VEP may return one value for a point overlap or ampersand-joined values for
    multi-base overlaps. Use the maximum finite segment value — the variant's
    functional impact is driven by the highest-expressed overlapping region.
    """
    if pd.isna(value):
        return np.nan
    if isinstance(value, (int, float, np.integer, np.floating)):
        return float(value)
    values = []
    for item in str(value).replace(",", "&").split("&"):
        item = item.strip()
        if item in {"", ".", "NA", "NaN", "nan"}:
            continue
        try:
            values.append(float(item))
        except ValueError:
            continue
    return max(values) if values else np.nan


def vep_consq_interpret_per_row(row: Dict) -> Tuple[bool, bool]:
    '''
    Evaluate the functional impact of a variant to the transcript from 2 perspectives:
    1. Whether the variant is a LoF variant
    2. Whether the variant is a protein length changing variant

    Args:
        row: A dictionary containing variant information

    Returns:
        Tuple[bool, bool]: (is_lof, is_length_changing)
    '''
    consq = row.get('Consequence', None)
    loftee_result = row.get('LoF', None)
    nmd_escaping = "escaping" in str(row.get('NMD', ""))
    if not isinstance(consq, str):
        return False, False

    # LOFTEE emits HC for conventional high-confidence pLoF calls and OS for
    # predicted "other splice" loss-of-function calls in extended splice sites
    # or at newly created donor sites. Both are high-confidence effect calls.
    loftee_tokens = {
        token.strip().upper()
        for token in re.split(r"[&,;|]", str(loftee_result or ""))
        if token.strip()
    }
    loftee_high_confidence = bool({"HC", "OS"} & loftee_tokens)

    lof_criteria = {
        'stop_gained',
        'start_lost',
        'transcript_ablation',
        'frameshift_variant'
    }

    # FIX: Added 'stop_lost' - stop-loss variants cause protein extension (PM4, not PVS1)
    # Per ClinGen SVI guidelines, stop-loss variants should be assigned PM4 as they cause
    # protein length change (extension) rather than loss-of-function
    length_changing_criteria = {
        'inframe_insertion',
        'inframe_deletion',
        'feature_elongation',
        'feature_truncation',
        'stop_lost'  # Stop-loss causes protein extension → PM4
    }

    consequence_lof = any(c in consq for c in lof_criteria)
    is_lof = loftee_high_confidence or (
        consequence_lof
        and "LC" not in loftee_tokens
        and not nmd_escaping
    )

    # For pure splicing related consequences, we trust SpliceAI and SpliceVault more than VEP
    is_length_changing = is_lof or any(c in consq for c in length_changing_criteria) or consequence_lof or nmd_escaping

    return is_lof, is_length_changing


def vep_consq_interpret(df: pd.DataFrame, threads: int = 10) -> pd.DataFrame:
    '''
    Apply the interpretation function to each row of the dataframe in parallel

    Args:
        df: Input DataFrame containing variant annotations
        threads: Number of CPU threads to use

    Returns:
        DataFrame with added 'lof' and 'len_changing' columns
    '''
    # Convert DataFrame to list of dicts for multiprocessing
    records = df.to_dict('records')

    # Create arguments for parallel processing
    with mp.Pool(threads) as pool:
        results = pool.map(vep_consq_interpret_per_row, records)

    # Add results as new columns
    df['vep_consq_lof'], df['vep_consq_len_changing'] = zip(*results)
    logger.info(f"vep_consq_interpret applied, {df['vep_consq_lof'].sum()} variants are having the LoF criteria")
    logger.info(f"vep_consq_interpret applied, {df['vep_consq_len_changing'].sum()} variants are having the protein length changing criteria")
    logger.info(f"vep_consq_interpret applied, now the table looks like: \n{df[:5].to_string(index=False)}")

    return df


def get_variant_type(hgvsp):
    """
    Determine the type of protein variant from HGVS notation.
    Returns 'nonsense', 'frameshift', 'delins', 'deletion', 'insertion_type' (for insertions and duplications),
    'missense', or 'unknown'.
    """
    if not isinstance(hgvsp, str):
        return None

    # Extract the p. part from potential ENSP prefixes (e.g., ENSP00000439902.1:p.Asn3264LeufsTer12)
    if ":" in hgvsp:
        hgvsp = hgvsp.split(":")[-1]

    # Order patterns based on biological severity
    # Check nonsense/stop gain variants first (most severe)
    if re.search(r'p\.\w+\d+Ter', hgvsp) or re.search(r'p\.\w+\d+\*', hgvsp):  # p.Arg123Ter or p.Arg123*
        return "nonsense"
    # Check frameshift second (next most severe)
    elif re.search(r'p\.\w+\d+fs', hgvsp) or re.search(r'p\.\w+\d+[Ff]s[Tt]er\d+', hgvsp):  # p.Gly123fs or p.Asn3264LeufsTer12
        return "frameshift"
    # Check delins (length changing)
    elif re.search(r'p\.\w+\d+delins', hgvsp):  # p.Ala123delinsGly
        return "delins"
    # Check deletion
    elif re.search(r'p\.\w+\d+del', hgvsp):  # p.Ala123del
        return "deletion"
    # Group duplications and insertions together as they're functionally similar (both add amino acids)
    elif re.search(r'p\.\w+\d+dup', hgvsp) or re.search(r'p\.\w+\d+ins', hgvsp):  # p.Ala123dup or p.Ala123insGly
        return "insertion_type"
    # Check missense last (least severe structural change) - more specific pattern to avoid false matches
    elif re.search(r'p\.[A-Z][a-z]+\d+[A-Z][a-z]+$', hgvsp):  # p.Ala123Gly (single amino acid change)
        return "missense"
    else:
        return None


def truncate_fraction(df):
    '''
    Calculate the fraction of the truncated protein for frameshift and stop_gained consequences
    '''
    nonsense_vars = df["Consequence"].str.contains("stop_gained") | df["Consequence"].str.contains("frameshift")
    # inframe_indels = df["Consequence"].str.contains("inframe_deletion") | df["Consequence"].str.contains("inframe_insertion")

    # Fix the issue by ensuring we maintain the same index
    protein_pos_series = df["Protein_position"].astype(str)
    contains_dash = protein_pos_series.str.contains("-", na=False)

    # Extract the first part for positions with "-"
    split_result = protein_pos_series.str.split("-", expand=True)
    first_part = split_result.iloc[:, 0] if len(split_result.columns) > 0 else protein_pos_series

    # Use pandas .where() instead of np.where() to maintain index alignment
    protein_poses = first_part.where(contains_dash, protein_pos_series)

    var_span = df["Protein_position"].map(lambda x: abs(pd.to_numeric(x.split("/")[0].split("-")[1], errors="coerce") - pd.to_numeric(x.split("/")[0].split("-")[0], errors="coerce") + 1) if "-" in str(x) else 1)
    var_protein_pos = pd.to_numeric(protein_poses.str.split("/", expand=True).iloc[:, 0], errors='coerce')
    total_protein_len = pd.to_numeric(df["Protein_position"].astype(str).str.split("/", expand=True).iloc[:, -1], errors='coerce')
    truncate_frac = np.where(nonsense_vars, 1 - (var_protein_pos / total_protein_len), var_span / total_protein_len)
    return truncate_frac
