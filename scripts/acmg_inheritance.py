#!/usr/bin/env python3
"""Legacy gene-level inheritance summaries and compatibility helpers.

The six variant-aware ACMG consumers no longer receive the arrays produced
here. They read the mechanism hub's condition-linked ``variant_inheritance``
and ``variant_penetrance`` columns. This module remains for the hub's audit
summary API and for callers that still request a coarse gene-level estimate.

    parse_hpo_inheritance              HPO inheritance terms -> one mode
    identify_inheritance_mode_per_row  one variant's mode, using the gene's
                                       mean AlphaMissense score and the
                                       ClinGen dosage call
    identify_inheritance_mode          the same over a whole table, in parallel
    hpo_onset_modes                    compatibility bridge using the shared
                                       penetrance vocabulary

The shared HPO vocabulary is deliberately independent of ACMG criteria, so the
hub can import this audit API without creating a criteria-module cycle.
"""

import logging
import re
import multiprocessing as mp
import pandas as pd
import numpy as np
from typing import Tuple

from hpo_penetrance import (
    HPO_INCOMPLETE_PENETRANCE_EQUIVALENT_TERMS,
    normalize_penetrance_hpo_ids,
)


logger = logging.getLogger(__name__)


def hpo_onset_modes(hpo_string):
    """Return whether unaffected-carrier observations are interpretable.

    This is a compatibility bridge for criteria that still read the old
    gene-wide HPO string.  It uses the same condition-modifier vocabulary as
    the mechanism hub: delayed or gradual onset and variable expressivity make
    an apparently unaffected carrier non-informative.  High penetrance,
    sex-limited expression, and slow progression are deliberately outside that
    set.

    Args:
        hpo_string: A string containing one or more HPO IDs
                    separated by semicolons. Can be None or empty.

    Returns:
        bool: False when an incomplete-penetrance-equivalent HPO term is
              present, otherwise True.

    **Disclaimer:** This check provides a simplified signal based on HPO terms
    related to onset and course. Always interpret gnomAD frequency using
    established guidelines (e.g., ACMG/AMP) and considering the specific
    disease context (prevalence, inheritance, penetrance, overall severity).
    """

    if not isinstance(hpo_string, str) or not hpo_string:
        return True

    hpo_ids = {
        term.strip().upper()
        for term in hpo_string.split(";")
        if re.match(r"^HP:\d+$", term.strip(), flags=re.IGNORECASE)
    }
    return not bool(hpo_ids & HPO_INCOMPLETE_PENETRANCE_EQUIVALENT_TERMS)


def parse_hpo_inheritance(row_dict: dict) -> str:
    # Parse the HPO_gene_inheritance field and return the inheritance mode
    # The HPO_gene_inheritance field is a string with multiple inheritance modes separated by semicolons
    # These inheritance modes can correspond to 3 different pathogenic mechanisms: LoF, GoF, DN.
    if isinstance(row_dict.get('HPO_IDs', None), str):
        hpo_terms = row_dict['HPO_IDs'].split(";")
        # This remains condition-coarse until all criteria consume the hub's
        # variant-condition histories, but its vocabulary must not disagree
        # with the higher-resolution path during that migration.
        incomplete_penetrance = (
            normalize_penetrance_hpo_ids(hpo_terms) == "incomplete"
        )
    else:
        incomplete_penetrance = False


    if isinstance(row_dict.get('HPO_gene_inheritance', None), str):
        hpo_inheritances = row_dict['HPO_gene_inheritance'].split(";")
    else:
        return incomplete_penetrance

    non_monogenic_set = {"Digenic inheritance", "Oligogenic inheritance", "Polygenic inheritance"}  # In most cases, these indicate compound heterozygous variants
    non_mendelian_set = {"Non-Mendelian inheritance"}  # Includes epigenetic modifications
    dominant_set = {"Autosomal dominant inheritance", "Autosomal dominant inheritance with maternal imprinting", "X-linked dominant inheritance"}
    # Treat generic X-linked inheritance as recessive by default. Male chrX
    # hemizygosity is handled later by sex-aware allele-state normalization.
    recessive_set = {"Autosomal recessive inheritance", "X-linked recessive inheritance", "X-linked inheritance"}

    # HPO recessive
    hpo_recessive = any([ hpo in recessive_set for hpo in hpo_inheritances ])
    # HPO dominant
    hpo_dominant = any([ hpo in dominant_set for hpo in hpo_inheritances ])
    # HPO non monogenic
    hpo_non_monogenic = any([ hpo in non_monogenic_set for hpo in hpo_inheritances ])
    # HPO non mendelian
    hpo_non_mendelian = any([ hpo in non_mendelian_set for hpo in hpo_inheritances ])

    return {
            'hpo_recessive': hpo_recessive,
            'hpo_dominant': hpo_dominant,
            'hpo_non_monogenic': hpo_non_monogenic,
            'hpo_non_mendelian': hpo_non_mendelian,
            'incomplete_penetrance': incomplete_penetrance
            }


def identify_inheritance_mode_per_row(row_dict: dict, gene_mean_am_score: float, clingen_curate_score: int = None) -> Tuple[bool, bool, bool, bool, bool, bool]:
    # We need to use three fields of the table to determine the inheritance mode:
    # 1. Gene
    # 2. LOEUF
    # 3. HPO_IDs
    # 4. HPO_gene_inheritance (overrides the above two fields), HPO observed dominant inheritance can derive from GOF variants
    # 5. ClinGen curated dosage sensitivity, 3 means haploinsufficient, 30 or 40 means haplosufficient

    loeuf_score = float(row_dict.get('LOEUF', 0.6))
    loeuf_score = 1.0 if pd.isna(loeuf_score) else loeuf_score  # If LOEUF is NaN, we leave the decision to gene avg AM score
    haplo_insufficient = (loeuf_score <= 0.35) or (gene_mean_am_score >= 0.564)
    haplo_insufficient = haplo_insufficient and ((loeuf_score <= 0.7) or pd.isna(loeuf_score)) and ((gene_mean_am_score >= 0.5) or pd.isna(gene_mean_am_score))
    haplo_sufficient = not haplo_insufficient

    clingen_recessive = None
    clingen_dominant = None
    if clingen_curate_score:
        logger.debug(f"Using ClinGen curated dosage sensitivity to determine inheritance mode for {row_dict['Gene']}, the ClinGen record looks like this: \n{clingen_curate_score}\n")
        if clingen_curate_score == 3:
            clingen_recessive = False
            clingen_dominant = True
            haplo_insufficient = True
        elif clingen_curate_score == 30 or clingen_curate_score == 40:
            clingen_recessive = True
            haplo_insufficient = False
            haplo_sufficient = True
            # Cannot modify hpo_dominant because it might relates to GOF variants

    hpo_inheritance = parse_hpo_inheritance(row_dict)
    if isinstance(hpo_inheritance, bool):
        logger.debug(f"No HPO inheritance information for {row_dict['Gene']}, using LOEUF: {row_dict['LOEUF']} and AM score: {gene_mean_am_score} to determine inheritance mode. The haploinsufficiency is {haplo_insufficient} and haplosufficiency is {haplo_sufficient}")
        recessive = clingen_recessive if clingen_recessive is not None else haplo_sufficient
        dominant = clingen_dominant if clingen_dominant is not None else haplo_insufficient
        return recessive, dominant, False, False, haplo_insufficient, hpo_inheritance

    if clingen_recessive is not None:
        hpo_inheritance['hpo_recessive'] = clingen_recessive

    if clingen_dominant is not None:
        hpo_inheritance['hpo_dominant'] = clingen_dominant

    if hpo_inheritance['hpo_recessive']:
        haplo_insufficient = False
        haplo_sufficient=True

    return hpo_inheritance['hpo_recessive'], hpo_inheritance['hpo_dominant'], hpo_inheritance['hpo_non_monogenic'], hpo_inheritance['hpo_non_mendelian'], haplo_insufficient, hpo_inheritance['incomplete_penetrance']


def identify_inheritance_mode(df: pd.DataFrame,
                              gene_to_am_score_map: dict,
                              clingen_dosage_sensitivity: str,
                              threads: int = 10) -> Tuple[np.ndarray, np.ndarray]:
    """
    Identify inheritance mode for each variant in parallel.

    Args:
        df: DataFrame containing variant information
        gene_to_am_score_map: Dictionary mapping genes to AM scores
        threads: Number of CPU threads to use

    Returns:
        Tuple of boolean arrays (dominant_array, recessive_array)
    """

    # Convert DataFrame rows to dictionaries for picklable input
    shrink_df = df.loc[:, ["Feature", "Gene", "SYMBOL", "LOEUF", "HPO_IDs", "HPO_gene_inheritance"]].drop_duplicates()
    row_dicts = shrink_df.to_dict('records')

    clingen_dosage_df = pd.read_table(clingen_dosage_sensitivity, low_memory=False).dropna(subset=["#Gene Symbol", "Haploinsufficiency Score"])
    clingen_dosage_map = dict(zip(clingen_dosage_df['#Gene Symbol'], clingen_dosage_df['Haploinsufficiency Score'].astype(int)))

    # Prepare arguments for starmap
    args = [(row_dict, gene_to_am_score_map.get(row_dict['Gene'], np.nan), clingen_dosage_map.get(row_dict['SYMBOL'], np.nan)) for row_dict in row_dicts]

    # Process in parallel using dictionaries instead of namedtuples
    threads = min(threads, len(row_dicts), mp.cpu_count()-1)
    with mp.Pool(threads) as pool:
        results = pool.starmap(identify_inheritance_mode_per_row, args)

    # Unzip results into separate arrays
    recessive_array, dominant_array, non_monogenic_array, non_mendelian_array, haplo_insufficient_array, incomplete_penetrance_array = zip(*results)
    shrink_df['recessive'] = np.array(recessive_array)
    shrink_df['dominant'] = np.array(dominant_array)
    shrink_df['non_monogenic'] = np.array(non_monogenic_array)
    shrink_df['non_mendelian'] = np.array(non_mendelian_array)
    shrink_df['haplo_insufficient'] = np.array(haplo_insufficient_array)
    shrink_df['incomplete_penetrance'] = np.array(incomplete_penetrance_array)

    # Map the arrays back to the original DataFrame, we need to use merge, anchor on Feature and Gene
    merged_df = df.merge(shrink_df, on=["Feature", "Gene", "SYMBOL", "LOEUF", "HPO_IDs", "HPO_gene_inheritance"], how="left")
    assert merged_df.shape[0] == df.shape[0], f"The number of rows in the merged DataFrame {merged_df.shape[0]} is not equal to the number of rows in the original DataFrame {df.shape[0]}"
    return merged_df.loc[:, "recessive"].fillna(False).to_numpy(), \
           merged_df.loc[:, "dominant"].fillna(False).to_numpy(), \
           merged_df.loc[:, "non_monogenic"].fillna(False).to_numpy(), \
           merged_df.loc[:, "non_mendelian"].fillna(False).to_numpy(), \
           merged_df.loc[:, "haplo_insufficient"].fillna(False).to_numpy(), \
           merged_df.loc[:, "incomplete_penetrance"].fillna(False).to_numpy()
