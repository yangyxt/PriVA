#!/usr/bin/env python3
"""Shared foundation: how the gene's disease is inherited, and when it starts.

No ACMG criterion lives here, but almost every benign-supporting one depends
on it. The arrays this module produces -- recessive, dominant, non-monogenic,
non-Mendelian, haploinsufficient, incomplete-penetrance -- are threaded by the
orchestrator into PP1, BS1, BS2, BS4, BP2/PM3 and BP5.

    parse_hpo_inheritance              HPO inheritance terms -> one mode
    identify_inheritance_mode_per_row  one variant's mode, using the gene's
                                       mean AlphaMissense score and the
                                       ClinGen dosage call
    identify_inheritance_mode          the same over a whole table, in parallel
    hpo_onset_modes                    age of onset, which BS2 needs to know
                                       whether a healthy carrier is old enough
                                       to be informative

This module imports nothing from the other acmg_* modules. That is what lets
gene_mechanism_hub.py import identify_inheritance_mode_per_row at the top of
the file instead of deferring the import inside a method to dodge a cycle.
"""

import logging
import re
import multiprocessing as mp
import pandas as pd
import numpy as np
from typing import Tuple


logger = logging.getLogger(__name__)


def hpo_onset_modes(hpo_string):
    """
    Checks if an HPO profile suggests early onset WITHOUT known confounding
    factors like late onset or slow/mild progression, which might affect
    gnomAD frequency interpretation.

    Args:
        hpo_string: A string containing one or more HPO IDs
                    separated by semicolons. Can be None or empty.

    Returns:
        bool: True if at least one 'early_onset' term is present AND
              no 'late_onset' or 'slow_mild' terms are present.
              False otherwise.

    **Disclaimer:** This check provides a simplified signal based on HPO terms
    related to onset and course. Always interpret gnomAD frequency using
    established guidelines (e.g., ACMG/AMP) and considering the specific
    disease context (prevalence, inheritance, penetrance, overall severity).
    """

    # --- Define HPO Sets ---

    # 1. List of Early Onset HPO IDs (Onset definitively before age 16)
    early_onset_hpos = {
        "HP:0030674",  # Antenatal onset (before birth)
        "HP:0011460",  # Embryonal onset (first 8 weeks)
        "HP:0011461",  # Fetal onset (after 8 weeks, before birth)
        "HP:0003577",  # Congenital onset (present at birth)
        "HP:0003623",  # Neonatal onset (<= 28 days)
        "HP:0003593",  # Infantile onset (28 days to 1 year)
        "HP:0011463",  # Childhood onset (1 to 5 years)
        "HP:0003621",  # Juvenile onset (5 to 15 years)
        "HP:0410280",  # Pediatric onset (broader term, 28 days to 15 years)
    }

    # 2. List of Late Onset HPO IDs (Onset at age 16 or later)
    late_onset_hpos = {
        "HP:0003596",  # Middle age onset (40-60 years)
        "HP:0003584",  # Late onset (>= 60 years)
    }
    # "HP:0003581",  # Adult onset (>= 16 years)
    # "HP:0011462",  # Young adult onset (16-40 years)

    # 3. List of Slow Progression or Mildness HPO IDs
    #    (Terms suggesting a course potentially compatible with survival/reproduction or reduced severity)
    slow_mild_hpos = {
        "HP:0031785",  # Insidious onset (gradual development)
        "HP:0012829",  # Mild (severity modifier)
        "HP:0040007",  # Asymptomatic
    }

    # Backup HPOs
    # "HP:0003774",  # Slow progression
    # "HP:0003678",  # Nonprogressive (condition does not worsen over time)
    # "HP:0011010",  # Chronic (persisting for a long time)

    # Combine the "confounding" factor lists for easier checking
    # These are terms that, if present, suggest caution is needed with gnomAD filtering
    confounding_hpos = late_onset_hpos.union(slow_mild_hpos)

    # Initialize flags
    found_early = False
    found_confounding = False

    if not isinstance(hpo_string, str) or not hpo_string:
        return True # Invalid input cannot meet criteria

    # Process input string
    potential_hpos = [term.strip() for term in hpo_string.split(';')]

    for hpo_id in potential_hpos:
        # Basic format check
        if re.match(r'^HP:\d+$', hpo_id):
            # Check if it's an early onset term
            if hpo_id in early_onset_hpos:
                found_early = True
            # Check if it's a confounding term (late onset OR slow/mild)
            if hpo_id in confounding_hpos:
                found_confounding = True
                # Optimization: If a confounding term is found, the final result
                # must be False, so we can stop iterating early.
                break

    # Evaluate the final condition:
    # Return True only if an early term was found AND no confounding term was found.
    return not found_confounding


def parse_hpo_inheritance(row_dict: dict) -> str:
    # Parse the HPO_gene_inheritance field and return the inheritance mode
    # The HPO_gene_inheritance field is a string with multiple inheritance modes separated by semicolons
    # These inheritance modes can correspond to 3 different pathogenic mechanisms: LoF, GoF, DN.
    if isinstance(row_dict.get('HPO_IDs', None), str):
        hpo_terms = row_dict['HPO_IDs'].split(";")
        # HP:0003829: Incomplete penetrance
        # HP:4000159: Moderate penetrance
        # HP:4000160: Low penetrance
        incomplete_penetrance = ("HP:0003829" in hpo_terms) or ("HP:4000159" in hpo_terms) or ("HP:4000160" in hpo_terms)
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
