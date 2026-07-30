#!/usr/bin/env python3
"""PS1 and PM5 -- this change, or this residue, has been classified before.

    PS1  the same amino acid change as a variant already established as
         pathogenic, reached by a different nucleotide change
    PM5  a different amino acid change at the same residue as an established
         pathogenic variant

Both are lookups against ClinVar at two stars or better. The vectorized path
matches on (transcript, protein position, HGVSp); check_splice_pathogenic is
the per-row fallback for variants with no protein consequence, which matches
on exon or intron position and SpliceAI agreement instead.
"""

import logging
import multiprocessing as mp
import pandas as pd
import numpy as np
from typing import Tuple

from splicing_var_analysis import parse_hgvsc_splice_position

from acmg_consequence import get_variant_type


logger = logging.getLogger(__name__)


def check_splice_pathogenic(row: dict,
                            clinvar_tranx_splice_dict_list: list,
                            pvs1_strength: int) -> bool:
    '''
    Check if a variant's splice change matches a known pathogenic variant.
    The dict is prepared by the script stat_aachange_clinvar.py
    The dict should actually be a list of dicts:
    pos_info = {
                'chrom': record.chrom,
                'pos': record.pos,
                'ref': record.ref,
                'alt': record.alts[0] if record.alts else None,
                'exon': exon,
                'intron': intron,
                'hgvsc': hgvsc,
                'consequence': consq,
                'clinvar_sig': cln_sig,
                'clinvar_review': rev_status,
                'splice_ai': splice_ai_data
            }
    while splice_ai_data is also a dict that looks like:
        SpliceAI_pred_DP_AG: <value> # Delta position for acceptor gain
        SpliceAI_pred_DP_AL: <value> # Delta position for acceptor loss
        SpliceAI_pred_DP_DG: <value> # Delta position for donor gain
        SpliceAI_pred_DP_DL: <value> # Delta position for donor loss
        SpliceAI_pred_DS_AG: <value> # Delta score for acceptor gain
        SpliceAI_pred_DS_AL: <value> # Delta score for acceptor loss
        SpliceAI_pred_DS_DG: <value> # Delta score for donor gain
        SpliceAI_pred_DS_DL: <value> # Delta score for donor loss
    '''
    if not clinvar_tranx_splice_dict_list:
        return False

    # Find the same splicing site pathogenic variant in the same transcript
    var_intron = row.get('INTRON', None)
    var_exon = row.get('EXON', None)
    if var_intron is None or var_exon is None:
        return False

    # Find the same splicing site pathogenic variant in the same transcript
    patho_records = [patho for patho in clinvar_tranx_splice_dict_list if patho['exon'] == var_exon or patho['intron'] == var_intron]

    if not patho_records:
        return False

    var_spliceai_ds_ag = row.get('SpliceAI_pred_DS_AG', np.nan)
    var_spliceai_ds_al = row.get('SpliceAI_pred_DS_AL', np.nan)
    var_spliceai_ds_dg = row.get('SpliceAI_pred_DS_DG', np.nan)
    var_spliceai_ds_dl = row.get('SpliceAI_pred_DS_DL', np.nan)

    delta_scores = abs(var_spliceai_ds_ag), abs(var_spliceai_ds_al), abs(var_spliceai_ds_dg), abs(var_spliceai_ds_dl)
    if all(ds < 0.5 for ds in delta_scores):
        return False

    # Find the biggest delta score index in delta_scores
    max_delta_score_index = np.argmax(delta_scores)
    max_delta_score = delta_scores[max_delta_score_index]
    target_spliceai_ds = ['SpliceAI_pred_DS_AG', 'SpliceAI_pred_DS_AL', 'SpliceAI_pred_DS_DG', 'SpliceAI_pred_DS_DL'][max_delta_score_index]
    logger.debug(f"The biggest DELTA score from SpliceAI is {max_delta_score}, which is from {target_spliceai_ds}")

    for patho in patho_records:
        patho_spliceai_ds = float(patho["splice_ai"].get(target_spliceai_ds, np.nan))
        if abs(patho_spliceai_ds) <= max_delta_score:
            vua_hgvsc = str(row["HGVSc"])
            vua_length = abs(len(row["ref"]) - len(row["alt"]))
            vua_length = 1 if vua_length == 0 else vua_length
            vua_parse_result = parse_hgvsc_splice_position(vua_hgvsc, str(row["STRAND"]), int(row["pos"]), vua_length)
            if vua_parse_result is None:
                vua_parse_result = {"overlapping_canonical_site": False}
            patho_length = abs(len(patho["ref"]) - len(patho["alt"]))
            patho_length = 1 if patho_length == 0 else patho_length
            patho_parse_result = parse_hgvsc_splice_position(patho["hgvsc"], str(row["STRAND"]), int(patho["pos"]), patho_length)
            if patho_parse_result is None:
                patho_parse_result = {"overlapping_canonical_site": False}
            patho_lp = "ikely" in patho["clinvar_sig"]

            if vua_parse_result["overlapping_canonical_site"]:
                logger.info(f"The variant {row['HGVSc']} is overlapping with a canonical splice site, the pathogenic variant is {patho['hgvsc']}")
                if pvs1_strength >= 4:
                    if patho_parse_result["overlapping_canonical_site"]:
                        logger.info(f"The pathogenic variant {patho['hgvsc']} is overlapping with a canonical splice site")
                        return False
                    else:
                        logger.info(f"The pathogenic variant {patho['hgvsc']} is not overlapping with a canonical splice site")
                        return "PS1_Supporting"
                elif pvs1_strength >= 1:
                    if patho_parse_result["overlapping_canonical_site"]:
                        logger.info(f"The pathogenic variant {patho['hgvsc']} is overlapping with a canonical splice site")
                        if patho_lp:
                            return False
                        else:
                            return "PS1"
                    else:
                        logger.info(f"The pathogenic variant {patho['hgvsc']} is not overlapping with a canonical splice site")
                        if patho_lp:
                            return "PS1_Supporting"
                        else:
                            return "PS1_Moderate"
                else:
                    return False
            else:
                logger.info(f"The variant {row['HGVSc']} is not overlapping with a canonical splice site, the pathogenic variant is {patho['hgvsc']}")
                if int(patho["pos"]) == int(row["pos"]) and str(row["chrom"]) == str(patho["chrom"]):
                    if patho_lp:
                        return "PS1_Moderate"
                    else:
                        return "PS1"
                else:
                    if patho_lp:
                        return "PS1_Supporting"
                    else:
                        return "PS1_Moderate"
    return False


def PS1_PM5_criteria(df: pd.DataFrame,
                     clinvar_aa_dict: dict,
                     clinvar_splice_dict: dict,
                     ps3_clinvar_patho: np.ndarray,
                     pvs1_criteria: np.ndarray,
                     threads: int = 10,
                     high_confidence_status = {
                                                'practice_guideline': 4,                                   # 4 stars
                                                'reviewed_by_expert_panel': 3,                             # 3 stars
                                                'criteria_provided,_multiple_submitters,_no_conflicts': 2,  # 2 stars
                                              }) -> Tuple[np.ndarray, np.ndarray]:
    '''
    Identify PS1/PM5 using vectorized pandas lookups against pre-flattened ClinVar sets,
    falling back to per-row check_splice_pathogenic only for rows without HGVSp.

    PS1: Same amino acid change as a previously established pathogenic variant.
    PM5: Same AA residue with a different missense change (same variant_type) as a previously
         established pathogenic variant.

    Both ClinVar structures arrive already loaded. The caller loads them once and
    shares them with summarize_clinvar_gene_pathogenicity, so the 44 MiB and
    211 MiB pickles are read from disk once per run rather than twice.
    '''
    # ---------------------------------------------------------------
    # Build flat lookup sets (one-time O(ClinVar_size) pass)
    #   ps1_set      = {(transcript, raw_protein_pos, hgvsp_alt)} for high-conf pathogenic
    #   pm5_residues = {(transcript, raw_protein_pos, variant_type)} for the same
    # We carry variant_type in the PM5 key because PM5 requires same variant_type
    # (missense vs nonsense vs frameshift) per the original per-row logic.
    # ---------------------------------------------------------------
    ps1_set = set()
    pm5_residues = set()
    for transcript, positions in clinvar_aa_dict.items():
        for raw_pos, variants in positions.items():
            for hgvsp_alt, entry in variants.items():
                is_patho_highconf = False
                for sig, rev_stat in zip(entry['CLNSIG'], entry['CLNREVSTAT']):
                    stars = high_confidence_status.get(rev_stat, 0)
                    if ("Pathogenic" in sig and stars >= 2) or \
                       (sig == "Likely_pathogenic" and stars >= 3):
                        is_patho_highconf = True
                        break
                if not is_patho_highconf:
                    continue
                ps1_set.add((transcript, raw_pos, hgvsp_alt))
                # Encode variant_type at the residue level so PM5 only matches same variant kind
                vtype = get_variant_type(hgvsp_alt)
                if vtype is not None:
                    pm5_residues.add((transcript, raw_pos, vtype))

    logger.info(f"Built PS1 lookup set with {len(ps1_set)} (transcript, pos, hgvsp) entries "
                f"and PM5 residue set with {len(pm5_residues)} (transcript, pos, variant_type) entries")

    # ---------------------------------------------------------------
    # Vectorized PS1: exact (Feature, Protein_position, HGVSp) tuple match
    # ---------------------------------------------------------------
    feature = df['Feature'].astype(str).values
    raw_pp = df['Protein_position'].astype(str).values
    hgvsp_col = df['HGVSp'].astype(str).values
    consequence = df['Consequence'].astype(str).values

    ps1_keys = list(zip(feature, raw_pp, hgvsp_col))
    ps1_eligible = np.array([k in ps1_set for k in ps1_keys], dtype=bool)

    # ---------------------------------------------------------------
    # Vectorized PM5: same (Feature, Protein_position, variant_type) but different HGVSp
    # PM5 only applies when variant_type == 'missense' for both query and ClinVar entry
    # (the original per-row logic excludes synonymous and requires same variant_type)
    # ---------------------------------------------------------------
    query_vtypes = np.array([get_variant_type(h) if h not in ('nan', 'inf', '') else None
                             for h in hgvsp_col], dtype=object)
    pm5_keys = list(zip(feature, raw_pp, query_vtypes))
    pm5_residue_match = np.array([k in pm5_residues for k in pm5_keys], dtype=bool)
    # PM5 requires same residue + same variant_type, but NOT same HGVSp (that would be PS1)
    pm5_eligible = pm5_residue_match
    # Exclude synonymous variants from PM5 (matches original logic at line 1034)
    synonymous_mask = pd.Series(consequence).str.contains('synonymous').fillna(False).values
    pm5_eligible = pm5_eligible & ~synonymous_mask & ~ps1_eligible

    logger.info(f"Vectorized lookup: PS1 hits={ps1_eligible.sum()}, PM5 hits={pm5_eligible.sum()}")

    # ---------------------------------------------------------------
    # Splice fallback: per-row logic for variants without HGVSp or where the
    # position is not in clinvar_aa_dict. The original per-row function returns
    # the splice result directly for these cases. Run only on the subset.
    # ---------------------------------------------------------------
    no_hgvsp = np.isin(hgvsp_col, ['nan', 'inf', ''])
    has_transcript_in_clinvar = np.array([f in clinvar_aa_dict or f in clinvar_splice_dict for f in feature], dtype=bool)
    pos_in_clinvar = np.array([(f in clinvar_aa_dict and p in clinvar_aa_dict[f])
                               for f, p in zip(feature, raw_pp)], dtype=bool)
    splice_fallback_mask = has_transcript_in_clinvar & (no_hgvsp | ~pos_in_clinvar)
    n_splice_fallback = splice_fallback_mask.sum()
    logger.info(f"Splice fallback subset: {n_splice_fallback} rows (no HGVSp or position absent in ClinVar AA dict)")

    splice_results = np.array([False] * len(df), dtype=object)
    if n_splice_fallback > 0:
        fallback_idx = np.where(splice_fallback_mask)[0]
        records_subset = df.iloc[fallback_idx].to_dict('records')
        args = [(records_subset[i],
                 clinvar_splice_dict.get(records_subset[i].get('Feature', ''), []),
                 pvs1_criteria[fallback_idx[i]]) for i in range(len(records_subset))]
        chunk_size = max(len(args) // (threads * 4), 1)
        with mp.Pool(threads) as pool:
            fallback_results = pool.starmap(check_splice_pathogenic, args, chunksize=chunk_size)
        for global_i, res in zip(fallback_idx, fallback_results):
            splice_results[global_i] = res

    # Splice fallback can return PS1-equivalent strength values
    splice_ps1_strong = (splice_results == "Same_AA_Change") | (splice_results == "PS1") | (splice_results == "PS1_PM5")
    splice_ps1_moderate = (splice_results == "PS1_Moderate")
    splice_ps1_supporting = (splice_results == "PS1_Supporting")
    splice_pm5 = (splice_results == "Same_AA_Residue") | (splice_results == "PS1_PM5")


    # Remove self-matching cases
    clinvar_lof = df['CLNSIG'].fillna("").str.contains('pathogenic') & (df['CLNREVSTAT'].map(high_confidence_status) >= 2) # Including Likely_pathogenic

    # ---------------------------------------------------------------
    # Combine vectorized + splice fallback results
    # ---------------------------------------------------------------
    ps1_array = np.zeros(len(df), dtype=int)
    pm5_array = np.zeros(len(df), dtype=int)

    # PS1 strong (3): vectorized hits OR splice fallback strong hits
    ps1_array[ps1_eligible | splice_ps1_strong] = 3
    # Splice-only moderate/supporting (only set if not already strong)
    ps1_array[splice_ps1_moderate & (ps1_array == 0)] = 2
    ps1_array[splice_ps1_supporting & (ps1_array == 0)] = 1
    ps1_array[clinvar_lof] = 0  # Remove self-matching for PS1 (LoF variants that are already pathogenic in ClinVar)

    pm5_eligible_combined = pm5_eligible | (splice_pm5 & ~synonymous_mask)
    # Remove PS1
    pm5_eligible_combined = pm5_eligible_combined & (ps1_array == 0)
    pm5_array[pm5_eligible_combined] = 2

    both_count = int(((ps1_array > 0) & (pm5_array > 0)).sum())
    logger.info(f"PS1/PM5 results (vectorized + splice fallback): "
                f"PS1 (any strength)={int((ps1_array > 0).sum())}, "
                f"PM5={int((pm5_array > 0).sum())}, "
                f"Both PS1+PM5={both_count}")
    logger.info(f"Self-matching prevention applied to {ps3_clinvar_patho.sum()} ClinVar pathogenic records")

    return ps1_array, pm5_array
