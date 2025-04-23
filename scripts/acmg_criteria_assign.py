#!/usr/bin/env python

import re
import pysam
import pickle
import logging
import pandas as pd
import numpy as np
import argparse as ap
from typing import Tuple, Dict
import multiprocessing as mp
from multiprocessing import Manager
from scipy.stats import binom
import mmap
import gc

from stat_protein_domain_amscores import nested_defaultdict
from protein_domain_mapping import DomainNormalizer
from mavdb_interpreter import MaveDBScoreInterpreter

logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)
console_handler=logging.StreamHandler()
console_handler.setLevel(logging.INFO)
formatter = logging.Formatter("%(levelname)s:%(asctime)s:%(funcName)s:%(lineno)s:%(message)s")
console_handler.setFormatter(formatter)
logger.addHandler(console_handler)


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
    if not isinstance(consq, str):
        return False, False
        
    # Skip splicing related consequences to leave them handled by SpliceAI and SpliceVault
    lof_criteria = {
        'stop_gained', 
        'stop_lost', 
        'start_lost', 
        'transcript_ablation', 
        'frameshift_variant'
    }
    
    length_changing_criteria = {
        'inframe_insertion',
        'inframe_deletion',
        'feature_elongation',
        'feature_truncation'
    }
    
    is_lof = any(c in consq for c in lof_criteria) | (loftee_result == 'HC')
    is_length_changing = is_lof or any(c in consq for c in length_changing_criteria)
    
    return is_lof, is_length_changing



def vep_consq_interpret(df: pd.DataFrame, threads: int = 10) -> pd.DataFrame:
    '''
    Apply the interpretation function to each row of the dataframe in parallel
    
    Args:
        df: Input DataFrame containing variant annotations
        threads: Number of CPU threads to use
        
    Returns:
        DataFrame with added 'lof' and 'length_changing' columns
    '''
    # Convert DataFrame to list of dicts for multiprocessing
    records = df.to_dict('records')
    
    # Create arguments for parallel processing
    with mp.Pool(threads) as pool:
        results = pool.map(vep_consq_interpret_per_row, records)
    
    # Add results as new columns
    df['vep_consq_lof'], df['vep_consq_length_changing'] = zip(*results)
    logger.info(f"vep_consq_interpret applied, {df['vep_consq_lof'].sum()} variants are having the LoF criteria")
    logger.info(f"vep_consq_interpret applied, {df['vep_consq_length_changing'].sum()} variants are having the protein length changing criteria")
    logger.info(f"vep_consq_interpret applied, now the table looks like: \n{df[:5].to_string(index=False)}")
    
    return df



def summarize_clinvar_gene_pathogenicity(clinvar_gene_aa_dict: dict) -> set:
    '''
    Based on the clinvar_gene_aa_dict:
    {ensg: {protein_pos: {hgvsp: {'CLNSIG': [cln_sig], 'CLNREVSTAT': [rev_stat]}}}}

    Extract the gene list that has been reported to harbor pathogenic variants in ClinVar database,
    considering both CLNSIG and CLNREVSTAT values.
    '''
    high_confidence_status = {
        'practice_guideline': 4,                                    # 4 stars
        'reviewed_by_expert_panel': 3,                              # 3 stars
        'criteria_provided,_multiple_submitters,_no_conflicts': 2,  # 2 stars
    }

    # Initialize set to store genes with high-confidence pathogenic variants
    pathogenic_genes = set()

    # Iterate through each gene and its amino acid changes
    for ensg, aa_positions in clinvar_gene_aa_dict.items():
        for pos, hgvsp_dict in aa_positions.items():
            for hgvsp, info in hgvsp_dict.items():
                # Retrieve CLNSIG and CLNREVSTAT lists
                clnsig_list = info.get('CLNSIG', [])
                revstat_list = info.get('CLNREVSTAT', [])

                # Check if any variant is pathogenic with high confidence
                for clnsig, revstat in zip(clnsig_list, revstat_list):
                    if ('Pathogenic' in clnsig) and (high_confidence_status.get(revstat, 0) == 2):
                        pathogenic_genes.add(ensg)
                        break  # No need to check other entries for this gene
                    if ('athogenic' in clnsig) and (high_confidence_status.get(revstat, 0) > 2):
                        pathogenic_genes.add(ensg)
                        break  # No need to check other entries for this gene
                if ensg in pathogenic_genes:
                    break  # Move to the next gene
            if ensg in pathogenic_genes:
                break  # Move to the next gene

    logger.info(f"Found {len(pathogenic_genes)} genes with high-confidence pathogenic variants in ClinVar")
    return pathogenic_genes


def identify_alternative_start_codon_genes(df: pd.DataFrame) -> set:
    if df["Consequence"].str.contains("start_lost").any():
        if df["Consequence"].str.contains("start_lost").all():
            return np.nan
        elif (np.logical_not(df["Consequence"].str.contains("start_lost")) & df["BIOTYPE"].str.contains("protein_coding")).any():
            return df["Gene"].values[0] if len(df["Gene"].values) > 0 else np.nan
        else:
            return np.nan
    else:
        return np.nan
    

def downstream_domain_impact(exon_str, tranx_id,tranx_exon_domain_map, interpro_entry_map_dict, dm_instance, domains=""):
    '''
    Used for frameshift and stopgain variants to explore whether the downstream protein region involving functional domains
    '''
    affected_exons = set([])
    if not isinstance(exon_str, str):
        pass
    elif "-" in exon_str and "/" in exon_str:
        affected_exons.update(range(int(exon_str.split("-")[0]), int(exon_str.split("/")[1]) + 1))
    elif "/" in exon_str:
        affected_exons.update(range(int(exon_str.split("/")[0]), int(exon_str.split("/")[1]) + 1))
    else:
        raise ValueError(f"Invalid exon string: {exon_str}")
    
    if isinstance(domains, str):
        domains = domains.split("&")
        for domain in domains:
            if dm_instance.interpret_functionality(domain, interpro_entry_map_dict) == "Functional":
                return True
    
    tranx_id = tranx_id.split(".", 1)[0] # Remove the ENSG version number
    for exon in affected_exons:
        if tranx_id in tranx_exon_domain_map:
            if exon in tranx_exon_domain_map[tranx_id]:
                domains = tranx_exon_domain_map[tranx_id][exon]
                for domain in domains:
                    domain = domain.split(":", 1)[1] # Remove the ENSG prefix in the domain path
                    if dm_instance.interpret_functionality(domain, interpro_entry_map_dict) == "Functional":
                        return True
    return False


def downstream_exon_patho_af(row, clinvar_patho_exon_af_dict, logic="any", threshold=0.01):
    exon_str = row['EXON']
    affected_exons = set([])
    if not isinstance(exon_str, str):
        return False
    elif "-" in exon_str and "/" in exon_str:
        affected_exons.update(range(int(exon_str.split("-")[0]), int(exon_str.split("/")[1]) + 1))
    elif "/" in exon_str:
        affected_exons.update(range(int(exon_str.split("/")[0]), int(exon_str.split("/")[1]) + 1))
    else:
        raise ValueError(f"Invalid exon string: {exon_str}")
    
    tranx_id = row['Feature']
    affected_exons_patho_af = set([])
    for exon in affected_exons:
        affected_exons_patho_af.add(clinvar_patho_exon_af_dict.get(tranx_id, {}).get(exon, (np.nan, ))[0])
    
    if logic == "any":
        return any(float(epa) < threshold for epa in affected_exons_patho_af)
    elif logic == "all":
        return all(float(epa) < threshold for epa in affected_exons_patho_af)
    else:
        raise ValueError(f"Invalid logic: {logic}, it should be either 'any' or 'all', depending on your needs")


def identify_functional_truncation(row, 
                                   dm_instance=None, 
                                   interpro_entry_map_dict=None, 
                                   tranx_exon_domain_map=None, 
                                   clinvar_patho_exon_af_dict=None, 
                                   exon_patho_af_threshold=0.01):
    '''
    Identify whether a truncating variant is involving functional regions on a protein
    '''
    if ("frameshift" in row['Consequence']) or ("stop_gained" in row['Consequence']) or row["5UTR_lof"] or row["splicing_lof"]:
        # These variants not only affect the local region of proteins, but also affect the downstream protein regions
        func_domain = downstream_domain_impact(row['EXON'], row['Feature'], tranx_exon_domain_map, interpro_entry_map_dict, dm_instance, domains=row["DOMAINS"])
        exon_frequent_patho = downstream_exon_patho_af(row, clinvar_patho_exon_af_dict, logic="any", threshold=exon_patho_af_threshold)
    else:
        if isinstance(row["DOMAINS"], str):
            func_domain = any(dm_instance.interpret_functionality(domain, interpro_entry_map_dict) == "Functional" for domain in row["DOMAINS"].split("&"))
        else:
            func_domain = False
        exon_frequent_patho = float(clinvar_patho_exon_af_dict.get(row['Feature'], {}).get(row['EXON'], (np.nan, ))[0]) < exon_patho_af_threshold
    
    return func_domain, exon_frequent_patho
    


def control_false_neg_rate(
    allele_frequencies: pd.Series,
    allele_numbers: pd.Series,
    af_threshold = 0.001,
    alpha: float = 0.01
) -> Tuple[pd.Series, pd.Series]:
    """
    Performs a one-sided binomial test for each variant with potentially
    different allele frequency thresholds for each variant.

    Args:
        allele_frequencies (pd.Series): Series of observed allele frequencies.
        allele_numbers (pd.Series): Series of total observed allele numbers.
        af_threshold (float, pd.Series, or np.ndarray): The allele frequency
            threshold(s) for the null hypothesis (H₀: p_true < af_threshold).
            If array-like, must match the length of allele_frequencies.
        alpha (float): The significance level for rejecting H₀. Default is 0.01.

    Returns:
        Tuple[pd.Series, pd.Series]: p_values and reject_h0 decisions.
    """
    if not isinstance(allele_frequencies, pd.Series):
        allele_frequencies = pd.Series(allele_frequencies)
    if not isinstance(allele_numbers, pd.Series):
        allele_numbers = pd.Series(allele_numbers)

    if not len(allele_frequencies) == len(allele_numbers):
        raise ValueError("Input Series must have the same length.")

    # Handle af_threshold input
    if isinstance(af_threshold, (float, int)):
        # If scalar, convert to Series with same index as allele_frequencies
        af_thresh_series = pd.Series(af_threshold, index=allele_frequencies.index)
    elif isinstance(af_threshold, np.ndarray):
        # Convert numpy array to Series, aligning index
        if len(af_threshold) != len(allele_frequencies):
            raise ValueError("If af_threshold is an array, it must have the same length as allele_frequencies.")
        af_thresh_series = pd.Series(af_threshold, index=allele_frequencies.index)
    elif isinstance(af_threshold, pd.Series):
        af_thresh_series = af_threshold
        if len(af_thresh_series) != len(allele_frequencies):
            raise ValueError("If af_threshold is a Series, it must have the same length as allele_frequencies.")
    else:
        raise TypeError("af_threshold must be a float, int, pandas Series, or numpy array.")
    
    # Ensure threshold is numeric and within [0,1]
    af_thresh_series = pd.to_numeric(af_thresh_series, errors='coerce').clip(0, 1)
    if af_thresh_series.isnull().any():
        raise ValueError("Non-numeric values found in af_threshold Series/array after coercion.")

    # --- Input Cleaning and Preparation (same as before) ---
    an = pd.to_numeric(allele_numbers, errors='coerce')
    an = an.fillna(0).astype(int)
    an[an < 0] = 0 

    af = pd.to_numeric(allele_frequencies, errors='coerce')
    af = af.clip(0, 1)

    ac_float = np.where(np.isnan(af), 0, af * an)
    ac_observed = np.where(np.isnan(ac_float), 0, np.round(ac_float).astype(int))
    ac_observed = np.minimum(ac_observed, an)
    ac_observed[an == 0] = 0

    valid_mask = (an > 0) & (~af.isna()) & (~allele_numbers.isna())

    # --- P-value Calculation ---
    p_values = pd.Series(np.nan, index=allele_frequencies.index)
    
    # Handle edge case k=0: P(AC >= 0) is always 1
    p_values.loc[valid_mask & (ac_observed == 0)] = 1.0
    
    # Calculate for k > 0
    mask_k_pos = valid_mask & (ac_observed > 0)
    if mask_k_pos.any():
        # Get the threshold values for each valid variant
        thresholds_to_use = af_thresh_series[mask_k_pos].values
        
        p_values.loc[mask_k_pos] = binom.sf(
            k=ac_observed[mask_k_pos] - 1,
            n=an[mask_k_pos],
            p=thresholds_to_use  # Now using the array of thresholds
        )

    # --- Decision ---
    reject_h0 = p_values <= alpha
    logger.info(f"There are {reject_h0.sum()} variants having their p-value less than {alpha}, {reject_h0.isna().sum()} datapoints are NA")

    return p_values, reject_h0


def PVS1_criteria(df: pd.DataFrame, 
                  am_score_dict: dict,
                  clinvar_gene_aa_dict: dict,
                  clinvar_patho_exon_af_stat: str,
                  interpro_entry_map_pkl: str,
                  ped_df: pd.DataFrame = None,
                  fam_name: str = None,
                  tranx_exon_domain_map_pkl: str = None) -> pd.DataFrame:
    # LoF is high confidence if it is a LoF variant and the consequence is known
    # VEP LoF, splicing LoF, UTRAnnotator LoF, ClinVar Pathogenic
    high_confidence_status = {
        # Higher confidence (2+ stars)
        'practice_guideline',                                    # 4 stars
        'reviewed_by_expert_panel',                             # 3 stars
        'criteria_provided,_multiple_submitters,_no_conflicts',  # 2 stars
    }

    clinvar_lof = (df['CLNSIG'] == 'Pathogenic') & df['CLNREVSTAT'].isin(high_confidence_status)
    logger.info(f"{clinvar_lof.sum()} variants are having ClinVar pathogenic variants")
    logger.info(f"{df['vep_consq_lof'].sum()} variants are considered LoF by VEP")
    logger.info(f"{df['splicing_lof'].sum()} variants are considered LoF by splicing effect prediction")
    logger.info(f"{df['5UTR_lof'].sum()} variants are considered LoF by 5UTR Annotator")
    lof_criteria = df['vep_consq_lof'] | df['splicing_lof'] | df['5UTR_lof'] | clinvar_lof

    # Do not consider the start_loss variants where functional transcripts for this gene uses alternative start codons
    alt_start_genes = set(df.groupby(['chrom', 'pos', 'ref', 'alt', 'Gene']).apply(identify_alternative_start_codon_genes).unique().tolist())
    logger.info(f"These are the {len(alt_start_genes)} genes that have functional transcripts using alternative start codons: {alt_start_genes}")
    alt_start_losts = df["Consequence"].str.contains("start_lost") & df["Gene"].isin(alt_start_genes)
    logger.info(f"{alt_start_losts.sum()} variants are having start_lost consequences to transcripts with alternative start codons")
    
    # Load the necessary dict file
    clinvar_patho_exon_af_dict = pickle.load(open(clinvar_patho_exon_af_stat, 'rb'))
    interpro_entry_map_dict = pickle.load(open(interpro_entry_map_pkl, 'rb'))
    tranx_exon_domain_map = pickle.load(open(tranx_exon_domain_map_pkl, 'rb'))
    dm_instance = DomainNormalizer()

    functional_domains, exon_rare_patho_afs = zip(*df.apply(identify_functional_truncation, axis=1, dm_instance=dm_instance, 
                                                                                                    interpro_entry_map_dict=interpro_entry_map_dict, 
                                                                                                    tranx_exon_domain_map=tranx_exon_domain_map, 
                                                                                                    clinvar_patho_exon_af_dict=clinvar_patho_exon_af_dict, 
                                                                                                    exon_patho_af_threshold=0.01))
    # Convert tuples to numpy arrays or pandas Series before applying ~
    functional_domains = np.array(functional_domains)  # Convert to numpy array
    exon_rare_patho_afs = np.array(exon_rare_patho_afs)  # Convert to numpy array
    unknown_function_domains = ~functional_domains
    logger.info(f"{unknown_function_domains.sum()} variants are truncating out protein regions with unknown functionality")
    exon_median_af = ~exon_rare_patho_afs
    logger.info(f"{exon_median_af.sum()} variants located in exons with median pathogenic AF > 0.01")
    nmd_escaping = df['NMD'].str.contains("NMD_escaping")
    logger.info(f"{nmd_escaping.sum()} variants are having NMD escaping, while {(nmd_escaping & lof_criteria).sum()} among them are predicted LoF")

    # stop-gain or frameshift variants that located in regions with unknown function
    stop_gains = df['Consequence'].str.contains("stop_gained")
    start_losts = df['Consequence'].str.contains("start_lost")
    frameshifts = df['Consequence'].str.contains("frameshift")
    inframe_dels = df['Consequence'].str.contains("inframe_deletion")
    inframe_dels = inframe_dels & ~stop_gains & ~start_losts & ~df['splicing_lof']

    null_variants_uncritical_region = unknown_function_domains & (stop_gains | frameshifts) & exon_median_af & nmd_escaping
    logger.info(f"{null_variants_uncritical_region.sum()} variants causing stop_gained or frameshift effects but not causing NMD, and only truncate out uncritical regions from the protein product")
    splicing_uncritical_region = ((df['splicing_len_changing'] | nmd_escaping) & unknown_function_domains & exon_median_af)
    logger.info(f"{splicing_uncritical_region.sum()} variants are having splicing altering, but did not cause NMD, and only cut off unknown functionality regions from protein")
    inframe_dels_uncritical_region = ((inframe_dels | nmd_escaping) & unknown_function_domains & exon_median_af)
    logger.info(f"{inframe_dels_uncritical_region.sum()} variants are having inframe deletions, did not cause NMD, and only cut off unknown functionality regions from protein")
    not_lof = null_variants_uncritical_region | splicing_uncritical_region | inframe_dels_uncritical_region | alt_start_losts
    logger.info(f"{not_lof.sum()} variants are not LoF, and located in exons with median pathogenic AF < 0.01")

    # Refine VEP LoF criteria
    df['vep_consq_length_changing'] = df['vep_consq_length_changing'] | ((null_variants_uncritical_region | inframe_dels_uncritical_region | splicing_uncritical_region) & df['vep_consq_lof'])
    df['vep_consq_lof'] = df['vep_consq_lof'] & ~(null_variants_uncritical_region | inframe_dels_uncritical_region | splicing_uncritical_region)
    logger.info(f"After refining,{df['vep_consq_lof'].sum()} variants are considered LoF by VEP")

    lof_criteria = lof_criteria & ~not_lof
    logger.info(f"In total, {lof_criteria.sum()} variants are considered LoF")

    '''
    Determine whether the gene is known to be pathogenic due to LoF
    Evaluated by mean AM scores and LOEUF
    1. if mean_AM_score is high and LoF depletion is high, then it means the gene is both functional essential and haploinsufficient
    2. if they are both low, then you cant tell whether the gene has functional redundancy or being haplosufficient.
    3. But for genes that are known to be pathogenic with AR inheritance mode, it is likely the gene is having low LoF depletion and low mean AM score and it is largely due to the haplo-sufficiency.
    
    There are two types of combinations we can use to determine.
    1. No matter the variant dosage (genotype), if the variant is LoF, and mean_AM_score is high or LoF depletion is high, then the variant is given PVS1.
    2. If the variant is LoF and homozygous, and the gene has been reported to be pathogenic in ClinVar, or the gene is with high mean_AM_score and LoF depletion, then the variant is given PVS1.
    '''

    mean_am_score_essential = df['Gene'].map(am_score_dict).fillna(0) >= 0.58
    loeuf_insufficient = df['LOEUF'] <= 0.35
    logger.info(f"{mean_am_score_essential.sum()} variants are having high mean AM score")
    logger.info(f"{loeuf_insufficient.sum()} variants are having low LOEUF, meaning intolerant to LoFs")
    intolerant_lof = mean_am_score_essential | loeuf_insufficient
    logger.info(f"In total, {intolerant_lof.sum()} variants are intolerant to LoFs based on mean AM score or LOEUF")

    clinvar_pathogenic_genes = summarize_clinvar_gene_pathogenicity(clinvar_gene_aa_dict)
    clinvar_pathogenic = df['Gene'].isin(clinvar_pathogenic_genes)
    logger.info(f"{clinvar_pathogenic.sum()} variants are having ClinVar pathogenic variants")
    
    clinvar_pathogenic = clinvar_pathogenic
    logger.info(f"{clinvar_pathogenic.sum()} variants are having ClinVar pathogenic variants and located in exons with median pathogenic AF < 0.01")

    if not ped_df is None and not fam_name is None:
        proband = ped_df.loc[(ped_df['#FamilyID'] == fam_name) & (ped_df['Phenotype'].isin(["2", 2])), 'IndividualID'].values[0]
        homozygous = df[proband].str.count("0") == 0
        pvs1_criteria = ( lof_criteria & intolerant_lof ) | \
                        ( lof_criteria & homozygous & ( intolerant_lof | clinvar_pathogenic ))
    else:
        logger.warning(f"No pedigree file provided, nor does we have a column to store the genotypes of the variants")
        pvs1_criteria = lof_criteria & ( intolerant_lof | clinvar_pathogenic )

    return pvs1_criteria



def check_splice_pathogenic(row: dict, 
                            clinvar_tranx_splice_dict_list: list) -> bool:
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
    
    for patho in patho_records:
        patho_spliceai_ds = patho.get(target_spliceai_ds, np.nan)
        if abs(patho_spliceai_ds) <= max_delta_score:
            return "Same_Splice_Site"
    
    return False


def extract_protein_position(hgvs_notation):
    # Pattern matches: "p." followed by letters, followed by numbers (position), followed by more letters
    pattern = r'p\.([A-Za-z]+)(\d+)([A-Za-z]+|\*)'
    
    match = re.search(pattern, hgvs_notation)
    if match:
        position = match.group(2)  # Group 2 contains just the position
        return int(position)
    else:
        return None


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
    elif re.search(r'p\.\w+\d+[A-Z][a-z]+$', hgvsp):  # p.Ala123Gly (single amino acid change)
        return "missense"
    else:
        return None


def check_aa_pathogenic(row: dict, 
                        clinvar_tranx_aa_dict: dict, 
                        clinvar_tranx_splice_dict_list: list, 
                        high_confidence_status: set) -> bool:
    '''
    Check if a variant's amino acid change matches a known pathogenic variant.
    
    Args:
        row: Variant annotation row
        clinvar_aa_dict: Nested dictionary containing ClinVar amino acid changes
        high_confidence_status: Set of high confidence review statuses
        
    Returns:
        bool: True if matches a pathogenic variant with high confidence
    '''
    transcript = row.get('Feature', '') # e.g ENST00000438441
    raw_protein_pos = row.get('Protein_position', '') # e.g 117/340 (pos/total_size)
    # protein_pos = str(raw_protein_pos).split("/")[0] if raw_protein_pos and not raw_protein_pos in [np.nan, np.inf] else '' # e.g "117"
    hgvsp = row.get('HGVSp', '') # e.g ENSP00000349098.5:p.E117K
    
    logger.debug(f"The current row records a variant overlapping with transcript {transcript} at protein position {raw_protein_pos} with HGVSp {hgvsp}")
        
    # Check if this transcript has any ClinVar entries
    if not clinvar_tranx_aa_dict:
        logger.debug(f"Transcript {transcript} not in ClinVar's VEP annotation records")
        return False
    
    if hgvsp in [np.nan, np.inf, 'nan', 'inf', '']:
        return check_splice_pathogenic(row, clinvar_tranx_splice_dict_list)
        
    # Check if this position has any ClinVar entries
    if raw_protein_pos not in clinvar_tranx_aa_dict:
        logger.debug(f"Protein position {raw_protein_pos} not in ClinVar's VEP annotation records for transcript {transcript}")
        return check_splice_pathogenic(row, clinvar_tranx_splice_dict_list)
        
    # Get the clinical significance and review status
    clinvar_entry = clinvar_tranx_aa_dict[raw_protein_pos].get(hgvsp, None)
    
    # Check if any entry is pathogenic with high confidence
    if clinvar_entry:
        logger.debug(f"There is a clinvar entry for {hgvsp} in transcript {transcript} at protein position {raw_protein_pos}, which is {clinvar_entry}")
        for sig, rev_stat in zip(clinvar_entry['CLNSIG'], clinvar_entry['CLNREVSTAT']):
            if ('athogenic' in sig) and (rev_stat in high_confidence_status):
                logger.debug(f"Same_AA_Change: {hgvsp} is pathogenic with high confidence in ClinVar")
                return "Same_AA_Change"
    
    logger.debug(f"No clinvar entry for {hgvsp} in transcript {transcript}. But there are AA changes recorded in the same protein position {raw_protein_pos}")
    splice_pathogenic = check_splice_pathogenic(row, clinvar_tranx_splice_dict_list)
    if splice_pathogenic:
        logger.debug(f"Same_Splice_Site: {hgvsp} is pathogenic with high confidence in ClinVar")
        return splice_pathogenic
    
    for hgvsp_alt, clinvar_entry in clinvar_tranx_aa_dict[raw_protein_pos].items():
        logger.debug(f"There is one clinvar entry for transcript {transcript} at protein position {raw_protein_pos}, which is variant {hgvsp_alt} has these clinvar annotations: {clinvar_entry}")
        hgvs_pos = extract_protein_position(hgvsp_alt)
        if hgvs_pos is None:
            logger.warning(f"The hgvs_alt {hgvsp_alt} is not a valid protein position, it is ignored")
            continue

        if str(hgvs_pos) != raw_protein_pos.split("/")[0]:
            logger.warning(f"The protein position {raw_protein_pos} does not match the hgvs_alt {hgvsp_alt}, it is ignored")
            continue
        
        # Check if variants are of the same type
        query_variant_type = get_variant_type(hgvsp)
        clinvar_variant_type = get_variant_type(hgvsp_alt)
        
        if query_variant_type != clinvar_variant_type:
            logger.debug(f"Variant types don't match: {hgvsp} ({query_variant_type}) vs {hgvsp_alt} ({clinvar_variant_type})")
            continue
        
        for sig, rev_stat in zip(clinvar_entry['CLNSIG'], clinvar_entry['CLNREVSTAT']):
            if ('athogenic' in sig) and (rev_stat in high_confidence_status):
                logger.debug(f"Same_AA_Residue: {hgvsp} is pathogenic with high confidence in ClinVar and is of the same type as {hgvsp_alt}")
                return "Same_AA_Residue"
            
    return False


def PS1_PM5_criteria(df: pd.DataFrame, 
                     clinvar_aa_dict_pkl: str, 
                     clinvar_splice_dict_pkl: str,
                     threads: int = 10) -> Tuple[np.ndarray, np.ndarray]:
    '''
    Identify variants using starmap,
    PS1: Same amino acid change as a previously established pathogenic variant
    PM5: Different amino acid change but same AA residue as a previously established pathogenic variant in a family member
    '''
    high_confidence_status = {
        'practice_guideline': 4,
        'reviewed_by_expert_panel': 3,
        'criteria_provided,_multiple_submitters,_no_conflicts': 2
    }
    logger.info(f"Loading ClinVar AA change dict from {clinvar_aa_dict_pkl}")
    clinvar_aa_dict = pickle.load(open(clinvar_aa_dict_pkl, 'rb'))
    logger.info(f"Loading ClinVar splice dict from {clinvar_splice_dict_pkl}")
    clinvar_splice_dict = pickle.load(open(clinvar_splice_dict_pkl, 'rb'))
    
    # Convert DataFrame to list of dictionaries
    records = df.to_dict('records')
    
    # Create argument tuples for starmap
    args = [(record, clinvar_aa_dict.get(record['Feature'], {}), clinvar_splice_dict.get(record['Feature'], {}), high_confidence_status) for record in records]
    
     # Add chunking
    chunk_size = max(len(records) // (threads * 4), 1)
    
    with mp.Pool(threads) as pool:
        logger.info(f"Running check_aa_pathogenic in parallel with {threads} threads on {len(records)} records")
        results = pool.starmap(check_aa_pathogenic, args, chunksize=chunk_size)
    
    results = np.array(results)
    ps1_criteria = (results == "Same_AA_Change") | (results == "Same_Splice_Site")
    pm5_criteria = results == "Same_AA_Residue"
    return ps1_criteria, pm5_criteria



def determine_denovo_per_row(row: dict, ped_info: dict, ped_df: pd.DataFrame) -> bool:
    # Input row is a dictionary converted from a pandas Series
    # Determine whether the variant is a denovo mutation in the proband
    proband_info, father_info, mother_info, sib_info = ped_info
    proband, proband_pheno = proband_info
    father, father_pheno = father_info
    mother, mother_pheno = mother_info

    proband_sex = ped_df.loc[ped_df['IndividualID'] == proband, 'Sex'].values[0]
    # First determine when proband is Male.
    if proband_sex == "1" or proband_sex == 1:
        if row['chrom'] == "chrX":
            # Variant on chrX are hemizygous in males
            if "1" in row.get(mother, ''):
                return False
            else:
                return "PS2"
        elif row['chrom'] == "chrY":
            # Variant on chrY are hemizygous in males
            if "1" in row.get(father, ''):
                return False
            else:
                return "PS2"
        elif row['chrom'] == "chrM":
            # Variant on chrM are hemizygous in males
            if "1" in row.get(mother, ''):
                return False
            else:
                return "PS2"
        else:
            if "1" in row.get(father, '') or "1" in row.get(mother, ''):
                return False
            elif father and mother:
                return "PS2"
            elif row.get('gnomAD_joint_AF') in [0, np.nan]:
                return "PM6"
            else:
                return False
    else:
        if row['chrom'] == "chrX":
            if "1" in row.get(father, '') or "1" in row.get(mother, ''):
                return False
            elif father and mother:
                return "PS2"
            elif row.get('gnomAD_joint_AF') in [0, np.nan]:
                return "PM6"
            else:
                return False
        elif row['chrom'] == "chrY":
            if "1" in row.get(father, ''):
                return False
            else:
                return "PS2"
        elif row['chrom'] == "chrM":
            if "1" in row.get(mother, ''):
                return False
            else:
                return "PS2"
        else:
            if "1" in row.get(father, '') or "1" in row.get(mother, ''):
                return False
            elif father and mother:
                return "PS2"
            elif row.get('gnomAD_joint_AF') in [0, np.nan]:
                return "PM6"
            else:
                return False



def PS2_PM6_criteria(df: pd.DataFrame, 
                     ped_df: pd.DataFrame, 
                     fam_name: str, 
                     threads: int = 10) -> Tuple[np.ndarray, np.ndarray]:
    """
    PS2: confirmed denovo mutation in the proband
    PM6: Assumed denovo mutation in the proband (if only have one parent and the PAF is absent from gnomAD)
    """
    logger.info(f"Running determine_denovo using pandas apply method")
    proband_info, father_info, mother_info, sib_info = identify_fam_members(ped_df, fam_name)
    proband, _ = proband_info
    father, _ = father_info
    mother, _ = mother_info
    
    # Get subset of pedigree for this family
    ped_subset = ped_df.loc[ped_df['#FamilyID'] == fam_name, :]
    
    # Select only necessary columns to reduce memory usage
    # These columns are needed for the denovo determination
    necessary_cols = ['chrom', 'pos', 'ref', 'alt', 'gnomAD_joint_AF']
    
    # Add columns for family members if they exist
    if proband in df.columns:
        necessary_cols.append(proband)
    if father in df.columns:
        necessary_cols.append(father)
    if mother in df.columns:
        necessary_cols.append(mother)
    
    # Create a unique key for each variant
    df['variant_key'] = df['chrom'] + '_' + df['pos'].astype(str) + '_' + df['ref'] + '_' + df['alt']
    
    # Keep track of original index
    original_index = df.index.copy()
    
    # Select necessary columns and drop duplicates
    cols_to_use = [col for col in necessary_cols if col in df.columns]
    unique_df = df[cols_to_use + ['variant_key']].drop_duplicates(subset=['variant_key'])
    
    logger.info(f"Processing {len(unique_df)} unique variants (reduced from {len(df)} total variants)")
    
    # Create ped_info tuple
    ped_info = (proband_info, father_info, mother_info, sib_info)
    
    # Apply the function to each unique row
    unique_results = unique_df.apply(
        lambda row: determine_denovo_per_row(row.to_dict(), ped_info, ped_subset), 
        axis=1
    )
    
    # Create a mapping from variant keys to results
    result_map = dict(zip(unique_df['variant_key'], unique_results))
    
    # Map results back to original dataframe
    results = np.array([result_map.get(key, False) for key in df['variant_key']])
    
    # Clean up temporary column
    df.drop('variant_key', axis=1, inplace=True)
    
    # Convert results to final form
    ps2_criteria = results == "PS2"
    pm6_criteria = results == "PM6"
    
    logger.info(f"Found {ps2_criteria.sum()} variants with PS2 criteria (confirmed de novo)")
    logger.info(f"Found {pm6_criteria.sum()} variants with PM6 criteria (assumed de novo)")
    
    return ps2_criteria, pm6_criteria




def mavedb_interpretation_per_row(row: pd.Series, urn_func_dict=None, mavedb_interpreter=None) -> bool:
    '''
    Interpret the MaveDB scores and return the PS3 and BS3 criteria
    '''
    score_sets = str(row.get('MaveDB_score', '')).split('&')
    high_confs = str(row.get('MaveDB_high_conf', '')).split('&')
    pvalues = str(row.get('MaveDB_pvalue', '')).split('&')
    scores = str(row.get('MaveDB_score', '')).split('&')
    urn_sets = str(row.get('MaveDB_urn')).split('&')

    score_interpretation = []
    mavedb_ps3 = False
    mavedb_bs3 = False
    for i, urn in enumerate(urn_sets):
        if urn not in urn_func_dict:
            logger.warning(f"The URN {urn} is not in the MaveDB metadata, it is ignored")
            continue
        interpretation = urn_func_dict[urn]
        score_interpretation.append(interpretation)
        score = scores[i]
        high_conf = high_confs[i]
        pvalue = pvalues[i]

        try:
            score = float(score)
        except ValueError:
            continue

        try:
            pvalue = float(pvalue)
        except ValueError:
            pvalue = np.nan

        result = mavedb_interpreter.interpret_score(score, interpretation, pvalue, high_conf)
        if result['mavedb_ps3']:
            mavedb_ps3 = True
        if result['mavedb_bs3']:
            mavedb_bs3 = True

    if mavedb_ps3 and mavedb_bs3:
        logger.warning(f"The variant {row['chrom']}:{row['pos']}:{row['ref']}:{row['alt']} has both PS3 and BS3 criteria according to MaveDB, the score is {score}, the high confidence is {high_conf}, the pvalue is {pvalue}, the interpretation is {interpretation}, it is likely to be a false positive")
        mavedb_ps3 = False
        mavedb_bs3 = False
        
    row["MaveDB_PS3"] = mavedb_ps3
    row["MaveDB_BS3"] = mavedb_bs3
    row["MaveDB_score_interpretation"] = "&".join([str(x) for x in score_interpretation])
    return row


def mavedb_score_interpretation(df: pd.DataFrame, mavedb_metadata: pd.DataFrame) -> pd.DataFrame:
    '''
    Interpret the MaveDB scores and return the PS3 and BS3 criteria
    '''
    urn_func_dict = mavedb_metadata.set_index('URN')['Score_Interpretation'].to_dict()
    mavedb_interpreter = MaveDBScoreInterpreter()
    df = df.apply(mavedb_interpretation_per_row, axis=1, urn_func_dict=urn_func_dict, mavedb_interpreter=mavedb_interpreter)
    mavedb_ps3_recs = df.loc[df['MaveDB_PS3'], [ "MaveDB_pvalue", "MaveDB_score", "MaveDB_high_conf", "MaveDB_score_interpretation"]]
    logger.info(f"There are {df['MaveDB_PS3'].sum()} variants with MaveDB determined PS3 criteria, they are determined by these functional assays: \n{mavedb_ps3_recs[:10].to_string()}")
    mavedb_bs3_recs = df.loc[df['MaveDB_BS3'], [ "MaveDB_pvalue", "MaveDB_score", "MaveDB_high_conf", "MaveDB_score_interpretation"]]
    logger.info(f"There are {df['MaveDB_BS3'].sum()} variants with MaveDB determined BS3 criteria, they are determined by these functional assays: \n{mavedb_bs3_recs[:10].to_string()}")
    return df


def PS3_BS3_criteria(df: pd.DataFrame, pvs1_criteria: np.ndarray, mavedb_metadata_tsv: str = "") -> pd.DataFrame:
    # Basically rely on ClinVar annotations
    high_confidence_status = {
        # Higher confidence (2+ stars)
        'practice_guideline': 4,                                   # 4 stars
        'reviewed_by_expert_panel': 3,                             # 3 stars
        'criteria_provided,_multiple_submitters,_no_conflicts': 2,  # 2 stars
    }
    
    clinvar_lof = df['CLNSIG'].str.contains('Pathogenic') & (df['CLNREVSTAT'].map(high_confidence_status) == 2)
    clinvar_lof = clinvar_lof | (df['CLNSIG'].str.contains('athogenic') & (df['CLNREVSTAT'].map(high_confidence_status) >= 3)) # Including Likely_pathogenic
    high_conf_benign = df['CLNSIG'].str.contains('Benign') & (df['CLNREVSTAT'].map(high_confidence_status) == 2)
    high_conf_benign = high_conf_benign | (df['CLNSIG'].str.contains('enign') & (df['CLNREVSTAT'].map(high_confidence_status) >= 3)) # Including Likely_benign

    clinvar_lof = clinvar_lof & ~pvs1_criteria  # Prevent inflation of matching a high revstatus patho record from clinvar
    
    if mavedb_metadata_tsv:
        mavedb_metadata = pd.read_table(mavedb_metadata_tsv, low_memory=False)
        mavedb_metadata.drop_duplicates(subset=["URN"], inplace=True)
        logger.info(f"There are {len(mavedb_metadata)} unique URNs in the MaveDB metadata, which looks like: \n{mavedb_metadata.head().to_string(index=False)}")
        df = mavedb_score_interpretation(df, mavedb_metadata)
        ps3_criteria = clinvar_lof | df['MaveDB_PS3']
        bs3_criteria = high_conf_benign | df['MaveDB_BS3']
    else:
        ps3_criteria = clinvar_lof
        bs3_criteria = high_conf_benign

    return {
            'PS3': ps3_criteria,
            'BS3': bs3_criteria
            }, df


def locate_intolerant_domain(row: dict, intolerant_domains: set) -> bool:
    '''
    Check if the variant is located in a mutational hotspot or a well-established functional protein domain
    The intolerant domains are calculated by scripts: 
    1. scripts/am_pick_intolerant_domains.py
    2. scripts/clinvar_pick_intolerant_domains.py
    '''
    gene = row.get('Gene', '') if isinstance(row.get('Gene', None), str) else ''
    if isinstance(row.get('DOMAINS', None), str):
        domains = [gene + ":" + d for d in row.get('DOMAINS', None).split('&')]
    else:
        return False

    return any(domain in intolerant_domains for domain in domains)



def fit_beta_mixture(x: np.ndarray) -> Tuple[float, bool]:
    """
    Fit a mixture of two beta distributions to determine if the distribution is bimodal.
    
    Args:
        x: Array of AM scores (between 0 and 1)
    Returns:
        Tuple[float, bool]: (Bimodality coefficient, Is_bimodal)
    """
    try:
        # Calculate basic statistics
        mean = np.mean(x)
        var = np.var(x)
        skewness = np.mean((x - mean) ** 3) / var ** 1.5
        kurtosis = np.mean((x - mean) ** 4) / var ** 2
        
        # Calculate bimodality coefficient
        # b = (skewness^2 + 1) / kurtosis
        # b > 0.555 indicates bimodality (empirical threshold)
        bimodality_coef = (skewness ** 2 + 1) / kurtosis
        
        return bimodality_coef, bimodality_coef > 0.555
        
    except Exception as e:
        logger.warning(f"Error in beta mixture fitting: {str(e)}")
        return 0.0, False


def analyze_score_distribution(scores: np.ndarray) -> Tuple[bool, float]:
    """
    Analyze if the AM score distribution is unimodal with a peak in pathogenic range.
    Only checks for unimodality and peak position, regardless of distribution symmetry.
    
    Args:
        scores: Array of AM scores
    
    Returns:
        Tuple[bool, float]: (is_pathogenic_unimodal, peak_score)
    """
    if len(scores) < 5:  # Need minimum points for analysis
        return False, 0.0
    
    try:
        # Create histogram
        hist, bin_edges = np.histogram(scores, bins=20, range=(0,1), density=True)
        bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2
        
        # Check for bimodality first
        _, is_bimodal = fit_beta_mixture(scores)
        
        if is_bimodal:
            return False, 0.0
            
        # If unimodal, find the peak
        peak_idx = np.argmax(hist)
        peak_score = bin_centers[peak_idx]
        
        # Check if peak is in pathogenic range (0.58-1.0)
        is_pathogenic = 0.58 <= peak_score <= 1.0
        
        return is_pathogenic, peak_score
        
    except Exception as e:
        logger.warning(f"Error in distribution analysis: {str(e)}")
        return False, 0.0


def locate_less_char_region(row: dict, am_intolerant_motifs: dict) -> Tuple[bool, bool]:
    '''
    Check if the variant is located in a less well-characterized region but with high AM scores. 
    The region was previous extracted during the installation step of the pipeline, it uses weighted KDE to identify the segments along the protein chain which are likely to be intolerant to AA changes
    '''

    aa_code_dict = {"Gly": "G", "Ala": "A", "Val": "V", "Leu": "L", 
                    "Ile": "I", "Thr": "T", "Ser": "S", "Met": "M", 
                    "Cys": "C", "Pro": "P", "Phe": "F", "Tyr": "Y", 
                    "Trp": "W", "His": "H", "Lys": "K", "Arg": "R", 
                    "Asp": "D", "Glu": "E", "Asn": "N", "Gln": "Q"}
    
    if not am_intolerant_motifs:
        logger.debug(f"No AM scores for gene {row['Gene']}, the AM profile might only calculate the variants effect at {row['chrom']}:{row['pos']} on another overlapping transcript that does not belong to this gene.")
        return False, False
    
    raw_protein_pos = row.get('Protein_position', '') # e.g 117/340 (pos/total_size)
    protein_pos = str(raw_protein_pos).split("/")[0] if raw_protein_pos and not raw_protein_pos in [np.nan, np.inf] else ''

    # Deal with non-protein-altering variants
    if not protein_pos:
        logger.debug(f"No protein position for variant {row['chrom']}:{row['pos']}, cannot determine if it is located in a mutational hotspot or a well-established functional protein domain")
        return False, False
    
    # Deal with indels
    if len(row['alt']) != len(row['ref']):
        return protein_pos in [x[1:] for x in am_intolerant_motifs['max_score_regions']], protein_pos in [x[1:] for x in am_intolerant_motifs['min_score_regions']]
    
    hgvsp_aa = row["HGVSp"].split(":")[1].lstrip("p.") if isinstance(row["HGVSp"], str) else ''
    # Extract the characters that is digit in the hgvsp_aa
    protein_position = re.search(r'\d+', hgvsp_aa).group() if hgvsp_aa else ''
    # Extract the string before the digits as the reference amino acid
    ref_aa = re.split(r'\d+', hgvsp_aa)[0] if protein_position else ''
    # Deal with stop codons
    if ref_aa == "Ter":
        logger.warning(f"The variant {row['chrom']}:{row['pos']} is a stop codon, can only determine by protein position")
        return protein_position in [x[1:] for x in am_intolerant_motifs['max_score_regions']], protein_position in [x[1:] for x in am_intolerant_motifs['min_score_regions']]
    # Deal with empty HGVSp
    if not ref_aa:
        return False, False
    
    single_letter_aa_code = aa_code_dict.get(ref_aa, ref_aa)
    assert len(single_letter_aa_code) == 1, f"For variant at {row['chrom']}:{row['pos']}, the amino acid code for {protein_position}: {ref_aa} is not a single letter"
    combo = f"{single_letter_aa_code}{protein_position}"
    logger.debug(f"The variant {row['chrom']}:{row['pos']} is causing AA changes at {combo}")
    key_regex = re.compile(rf'^[A-Z]{protein_position}$')

    return (combo in am_intolerant_motifs['max_score_regions']) or any(key_regex.match(x) for x in am_intolerant_motifs['max_score_regions']), \
           (combo in am_intolerant_motifs['min_score_regions']) or any(key_regex.match(x) for x in am_intolerant_motifs['min_score_regions'])

    


def PM1_criteria(df: pd.DataFrame, 
                 intolerant_domains_pkl: str, 
                 intolerant_motifs_pkl: str,
                 threads: int = 10) -> np.ndarray:
    # Memory-mapped approach
    with open(intolerant_domains_pkl, 'rb') as f:
        mm = mmap.mmap(f.fileno(), 0, access=mmap.ACCESS_READ)
        intolerant_domains = pickle.load(mm)

    logger.info(f"Loading the recorded intolerant domains which look alike: {intolerant_domains}")
    
    row_dicts = df.to_dict('records')
    args = [(row, intolerant_domains) for row in row_dicts]
    
    with mp.Pool(threads) as pool:
        results = pool.starmap(locate_intolerant_domain, args)

    loc_intol_domain = np.array(results)
    logger.info(f"There are {np.sum(loc_intol_domain)} variants located in a protein domain that is seemingly intolerant to AA changes according to AM scores")
    
    truncating = df['vep_consq_length_changing'] | df['splicing_len_changing'] | df['5UTR_len_changing'] | df['vep_consq_lof'] | df['splicing_lof'] | df['5UTR_lof']

    missense = df['Consequence'].str.contains('missense_variant')
    logger.info(f"There are {missense.sum()} missense variants")
    missense_damaging = df["am_class"].str.contains('athogenic')
    logger.info(f"There are {missense_damaging.sum()} missense variants that are considered damaging by AlphaMissense")
    missense_damaging = missense_damaging | (df["PrimateAI"] > 0.9)
    logger.info(f"There are {missense_damaging.sum()} missense variants that are considered damaging after considering the PrimateAI score")

    missense_benign = df["am_class"].str.contains('benign')
    if "MaveDB_PS3" in df.columns:
        # If DMS says the variant is damaging (decreased activity or expression), then consider it as damaging
        missense_damaging = missense_damaging | df["MaveDB_PS3"]
        logger.info(f"There are {missense_damaging.sum()} missense variants that are considered damaging after considering the MaveDB DMS/MPRA data")
    missense_damaging = missense_damaging & missense

    intolerant_motifs = pickle.load(open(intolerant_motifs_pkl, 'rb'))
    args = [(row, intolerant_motifs.get(row['Feature'], {})) for row in row_dicts]
    logger.info(f"There are {len(args)} variants to be checked for intolerant motifs")
    with mp.Pool(threads) as pool:
        results = pool.starmap(locate_less_char_region, args)
        # Results are a list of tuples, we need to convert them to two independent boolean arrays
        max_am_score_motifs = np.array([r[0] for r in results])
        min_am_score_motifs = np.array([r[1] for r in results])

    logger.info(f"There are {np.sum(min_am_score_motifs)} variants located in a mutational hotspot that is seemingly intolerant to AA changes according to AM scores")
    logger.info(f"There are {np.sum(max_am_score_motifs)} variants located in a mutational hotspot that is seemingly intolerant to AA changes according to most severe AM scores")
    return ( max_am_score_motifs & (missense_damaging | truncating) ) | \
           ( min_am_score_motifs & (truncating | missense) & ~missense_benign ) | \
           ( loc_intol_domain & (missense | truncating) & ~missense_benign )



def PM2_criteria(df: pd.DataFrame, gnomAD_extreme_rare_threshold: float = 0.0001) -> np.ndarray:
    # PM2: The variant is absent from gnomAD or the variant is extremely rare in gnomAD
    gnomAD_absent = (df['gnomAD_joint_AF'] == 0) | (df['gnomAD_joint_AF'].isna())
    # _, gnomAD_max_not_that_rare = control_false_neg_rate(df['gnomAD_joint_AF_max'], df['gnomAD_joint_AN_max'], af_threshold=gnomAD_extreme_rare_threshold, alpha=0.01)
    gnomAD_rare = df['gnomAD_joint_AF'] < gnomAD_extreme_rare_threshold
    return gnomAD_absent | gnomAD_rare


def PM4_criteria(df: pd.DataFrame) -> np.ndarray:
    # PM4: The variant is causing the protein length change
    return df['vep_consq_length_changing'] | df['splicing_len_changing'] | df['5UTR_len_changing']



def PP2_BP1_criteria(df: pd.DataFrame, 
                     domain_mechanism_tsv: str, 
                     conf_level: float = 0.05,
                     chunk_size: int = 10000) -> Tuple[pd.Series, pd.Series]:
    """
    Memory-efficient implementation of PP2/BP1 criteria evaluation.
    
    Args:
        df: Variant annotation DataFrame
        domain_mechanism_tsv: Path to domain mechanism analysis results
        conf_level: P-value threshold for significance
        
    Returns:
        Tuple[pd.Series, pd.Series]: PP2 and BP1 criteria boolean series
    """
    # Read domain mechanism data more efficiently by selecting only needed columns
    domain_df = pd.read_table(
        domain_mechanism_tsv, 
        usecols=['gene_id', 'domain', 'pvalue', 'gene_pvalue'],
        low_memory=False
    )
    
    # Create sets of tolerant/intolerant genes and domains using boolean indexing
    # This avoids creating intermediate DataFrames
    tolerant_genes = set(domain_df.loc[domain_df['gene_pvalue'] > conf_level, 'gene_id'])
    intolerant_genes = set(domain_df.loc[domain_df['gene_pvalue'] < conf_level, 'gene_id'])
    tolerant_domains = set(domain_df.loc[domain_df['pvalue'] > conf_level, 'domain'])
    intolerant_domains = set(domain_df.loc[domain_df['pvalue'] < conf_level, 'domain'])

    logger.info(f"There are {len(tolerant_genes)} tolerant genes and {len(intolerant_genes)} intolerant genes")
    logger.info(f"There are {len(tolerant_domains)} tolerant domains and {len(intolerant_domains)} intolerant domains")
    
    # Clear domain_df from memory
    del domain_df
    gc.collect()
    
    # Initialize result arrays
    n_variants = len(df)
    pp2_criteria = np.zeros(n_variants, dtype=bool)
    bp1_criteria = np.zeros(n_variants, dtype=bool)
    
    # Process in chunks to reduce memory usage
    for start_idx in range(0, n_variants, chunk_size):
        end_idx = min(start_idx + chunk_size, n_variants)
        chunk_slice = slice(start_idx, end_idx)
        
        # Process only variants that are missense
        missense_mask = df.iloc[chunk_slice]['Consequence'].str.contains('missense_variant')
        if not missense_mask.any():
            continue
            
        # Get relevant data for missense variants
        chunk_genes = df.iloc[chunk_slice]['Gene'][missense_mask]
        chunk_domains = df.iloc[chunk_slice]['DOMAINS'][missense_mask]
        
        # Process domains for missense variants
        for i, (gene, domains_str) in enumerate(zip(chunk_genes, chunk_domains)):
            if pd.isna(domains_str):
                continue
                
            # Create domain identifiers
            domain_ids = [f"{gene}:{domain}" for domain in str(domains_str).split('&')]
            
            # Check domain matches
            has_tolerant_domain = any(d in tolerant_domains for d in domain_ids)
            has_intolerant_domain = any(d in intolerant_domains for d in domain_ids)
            
            # Update criteria arrays for this variant
            true_idx = chunk_slice.start + missense_mask[:i+1].sum() - 1
            pp2_criteria[true_idx] = (gene in tolerant_genes) or has_tolerant_domain
            bp1_criteria[true_idx] = (gene in intolerant_genes) or has_intolerant_domain
    
    return pd.Series(pp2_criteria, index=df.index), pd.Series(bp1_criteria, index=df.index)



def PP3_BP4_criteria(df: pd.DataFrame) -> pd.Series:
    # PP3: predicted to be deleterious by in-silico tools
    # Including PrimateAI (Missense), CADD, AlphaMissense (Missense), VEP, SpliceAI, SpliceVault, UTRAnnotator

    # BP4: variant is reported benign
    return (df['PrimateAI'] > 0.9) | \
           (df['CADD_phred'] >= 20) | \
           (df['am_class'].str.contains('pathogenic')) | \
           (df['vep_consq_lof'] == True) | \
           (df['splicing_lof'] == True) | \
           (df['5UTR_lof'] == True), \
           (df['PrimateAI'] < 0.8) & \
           (df['CADD_phred'] < 20) & \
           (np.logical_not(df['am_class'].str.contains('athogenic'))) & \
           (df['vep_consq_lof'] != True) & \
           (df['splicing_lof'] != True) & \
           (df['5UTR_lof'] != True)


def PP5_BP6_criteria(df: pd.DataFrame) -> pd.Series:
    # PP5: The variant is reported as pathogenic by a reputable source but without to many supporting evidences
    return df['CLNSIG'].str.contains('athogenic'), \
           df['CLNSIG'].str.contains('enign')


def BS1_criteria(df: pd.DataFrame, 
                 gene_to_am_score_map: dict,
                 expected_incidence: float = 0.001, 
                 clinvar_patho_af_stat: str = "",
                 threads: int = 10):
    # BS1: PAF of variant is greater than expected incidence of the disease
    autosomal = (df['chrom'] != "chrX") & (df['chrom'] != "chrY")
    x_linked = df['chrom'] == "chrX"
    y_linked = df['chrom'] == "chrY"

    recessive, dominant, non_monogenic, non_mendelian, haplo_insufficient, incomplete_penetrance = identify_inheritance_mode(df, gene_to_am_score_map, threads)

    false_neg_rate, common_vars = control_false_neg_rate(df['gnomAD_joint_AF_max'], df['gnomAD_joint_AN_max'], af_threshold=expected_incidence, alpha=0.01)

    max_ind_incidence = np.where(df['gnomAD_joint_AN_max']/2 > 10/expected_incidence, df['gnomAD_nhomalt_max']/(df['gnomAD_joint_AN_max']/2), (df['gnomAD_nhomalt_XX'] + df['gnomAD_nhomalt_XY'])/(df['gnomAD_joint_AN']/2))
    max_af_larger_incidence = np.where(common_vars.isna(), df['gnomAD_joint_AF'] > expected_incidence, common_vars)
    logger.info(f"There are {max_af_larger_incidence.sum()} variants having their PAF greater than the expected incidence of the disease")
    # For autosomal dominant disease, we can assign BS1 if the variant is observed the frequency of the variant is greater than the expected incidence of the disease
    autosomal_dominant = autosomal & dominant & max_af_larger_incidence & np.logical_not(recessive) & np.logical_not(non_monogenic) & np.logical_not(non_mendelian)
    autosomal_recessive = autosomal & recessive & (max_ind_incidence > expected_incidence) & np.logical_not(non_monogenic) & np.logical_not(non_mendelian)

    # For X-linked disease, we can assign BS1 if the variant is observed the frequency of the variant is greater than the expected incidence of the disease
    x_linked_recessive = x_linked & recessive & ((df['gnomAD_nhomalt_XX'] + df['gnomAD_nhomalt_XY'])/(df['gnomAD_joint_AN']/2) > expected_incidence) & np.logical_not(non_monogenic) & np.logical_not(non_mendelian)
    x_linked_dominant = x_linked & dominant & max_af_larger_incidence & np.logical_not(non_monogenic) & np.logical_not(non_mendelian)
    # For Y-linked disease, we can assign BS1 if the variant is observed the frequency of the variant is greater than the expected incidence of the disease
    y_linked = y_linked & max_af_larger_incidence & np.logical_not(non_monogenic) & np.logical_not(non_mendelian)
    greater_than_disease_incidence = autosomal_dominant | autosomal_recessive | x_linked_recessive | x_linked_dominant | y_linked
    logger.info(f"There are {greater_than_disease_incidence.sum()} variants having their PAF greater than the expected incidence of the disease")
    
    logger.info(f"Loading the clinvar pathogenic AF stat from {clinvar_patho_af_stat}")
    clinvar_patho_af_dict = pickle.load(open(clinvar_patho_af_stat, 'rb'))
    clinvar_patho_af_dict = {k: v for k, v in clinvar_patho_af_dict.items() if v is not None}
    logger.info(f"The clinvar pathogenic AF stat type is {type(clinvar_patho_af_dict)}")
    df.loc[:, "Gene"] = df["Gene"].fillna(np.nan)
    gene_max_patho_af = df['Gene'].map(lambda gene: clinvar_patho_af_dict.get(gene, {"af": 0}).get('af', 0))
    _, greater_than_clinvar_patho_af = control_false_neg_rate(df['gnomAD_joint_AF_max'], df['gnomAD_joint_AN_max'], af_threshold=gene_max_patho_af, alpha=0.01)
    greater_than_clinvar_patho_af = np.where(greater_than_clinvar_patho_af.isna(), df['gnomAD_joint_AF'] > gene_max_patho_af, greater_than_clinvar_patho_af)
    _, greater_than_basic_af = control_false_neg_rate(df['gnomAD_joint_AF_max'], df['gnomAD_joint_AN_max'], af_threshold=0.0001, alpha=0.01)
    greater_than_basic_af = np.where(greater_than_basic_af.isna(), df['gnomAD_joint_AF'] > 0.0001, greater_than_basic_af)
    return (greater_than_disease_incidence | greater_than_clinvar_patho_af) & greater_than_basic_af, recessive, dominant, non_monogenic, non_mendelian, haplo_insufficient, incomplete_penetrance


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
        "HP:0003581",  # Adult onset (>= 16 years)
        "HP:0011462",  # Young adult onset (16-40 years)
        "HP:0003596",  # Middle age onset (40-60 years)
        "HP:0003584",  # Late onset (>= 60 years)
    }

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
        return False # Invalid input cannot meet criteria

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
    return found_early and not found_confounding


def parse_hpo_inheritance(row_dict: dict) -> str:
    # Parse the HPO_gene_inheritance field and return the inheritance mode
    # The HPO_gene_inheritance field is a string with multiple inheritance modes separated by semicolons
    # These inheritance modes can correspond to 3 different pathogenic mechanisms: LoF, GoF, DN. 
    if isinstance(row_dict.get('HPO_IDs', None), str):
        hpo_terms = row_dict['HPO_IDs'].split(";")
        incomplete_penetrance = "HP:0003829" in hpo_terms
    else:
        incomplete_penetrance = False


    if isinstance(row_dict.get('HPO_gene_inheritance', None), str):
        hpo_inheritances = row_dict['HPO_gene_inheritance'].split(";")
    else:
        return incomplete_penetrance
    
    non_monogenic_set = {"Digenic inheritance", "Oligogenic inheritance", "Polygenic inheritance"}  # In most cases, these indicate compound heterozygous variants
    non_mendelian_set = {"Non-Mendelian inheritance"}  # Includes epigenetic modifications
    dominant_set = {"Autosomal dominant inheritance", "Autosomal dominant inheritance with maternal imprinting", "X-linked dominant inheritance"}
    recessive_set = {"Autosomal recessive inheritance", "X-linked recessive inheritance"}


    # HPO recessive
    hpo_recessive = any( hpo in recessive_set for hpo in hpo_inheritances )
    # HPO dominant
    hpo_dominant = any( hpo in dominant_set for hpo in hpo_inheritances )
    # HPO non monogenic
    hpo_non_monogenic = any( hpo in non_monogenic_set for hpo in hpo_inheritances )
    # HPO non mendelian
    hpo_non_mendelian = any( hpo in non_mendelian_set for hpo in hpo_inheritances )

    return {
            'hpo_recessive': hpo_recessive,
            'hpo_dominant': hpo_dominant,
            'hpo_non_monogenic': hpo_non_monogenic,
            'hpo_non_mendelian': hpo_non_mendelian,
            'incomplete_penetrance': incomplete_penetrance
            }


def identify_inheritance_mode_per_row(row_dict: dict, gene_mean_am_score: float) -> Tuple[bool, bool]:
    # We need to use three fields of the table to determine the inheritance mode:
    # 1. Chromosome
    # 2. LOEUF
    # 3. AM score
    # 4. HPO_gene_inheritance (overrides the above two fields), HPO observed dominant inheritance can derive from GOF variants

    haplo_insufficient = (float(row_dict.get('LOEUF', np.nan)) <= 0.6) | (gene_mean_am_score >= 0.58)
    haplo_sufficient = not haplo_insufficient

    hpo_inheritance = parse_hpo_inheritance(row_dict)
    if isinstance(hpo_inheritance, bool):
        logger.debug(f"No HPO inheritance information for {row_dict['Gene']}, using LOEUF and AM score to determine inheritance mode. The row looks like this: \n{row_dict}\n")
        return haplo_insufficient, haplo_sufficient, False, False, haplo_insufficient, hpo_inheritance

    return hpo_inheritance['hpo_recessive'], hpo_inheritance['hpo_dominant'], hpo_inheritance['hpo_non_monogenic'], hpo_inheritance['hpo_non_mendelian'], haplo_insufficient, hpo_inheritance['incomplete_penetrance']

    

def identify_inheritance_mode(df: pd.DataFrame, 
                              gene_to_am_score_map: dict, 
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
    row_dicts = df.to_dict('records')
    
    # Prepare arguments for starmap
    args = [(row_dict, gene_to_am_score_map.get(row_dict['Gene'], np.nan)) for row_dict in row_dicts]
    
    # Process in parallel using dictionaries instead of namedtuples
    threads = min(threads, len(row_dicts), mp.cpu_count()-1)
    with mp.Pool(threads) as pool:
        results = pool.starmap(identify_inheritance_mode_per_row, args)
    
    # Unzip results into separate arrays
    recessive_array, dominant_array, non_monogenic_array, non_mendelian_array, haplo_insufficient_array, incomplete_penetrance_array = zip(*results)
    return np.array(recessive_array), np.array(dominant_array), np.array(non_monogenic_array), np.array(non_mendelian_array), np.array(haplo_insufficient_array), np.array(incomplete_penetrance_array)



def BS2_criteria(df: pd.DataFrame, 
                 bs1_criteria: np.ndarray,
                 recessive: np.ndarray,
                 dominant: np.ndarray,
                 non_monogenic: np.ndarray,
                 non_mendelian: np.ndarray,
                 haplo_insufficient: np.ndarray,
                 incomplete_penetrance: np.ndarray) -> pd.Series:
    '''
    The BS2 criteria is about observing variant in healthy adult
    There are several categories of situations here:
    1. For gene that is haplo-insufficient (autosomal dominant), we can assign BS2 if the variant is observed in a healthy adult (either homozygous or heterozygous)
    2. For gene that is autosomal recessive, we can assign BS2 if the variant is observed in a healthy adult with a homozygous genotype
    3. For X-linked recessive disease (gene on X but is haplo-sufficient), we can assign BS2 if the variant is observed in a healthy adult male (hemizygous) or a healthy adult female (homozygous)
    4. For X-linked dominant disease (gene on X and is haplo-insufficient), we can assign BS2 if the variant is observed in a healthy adult male (hemizygous)

    For haplo-insufficiency, we need to use LOEUF and mean AM score to determine.
    For panetrance, we use hpo terms to determine: HP:0003829
    '''
    autosomal = (df['chrom'] != "chrX") & (df['chrom'] != "chrY")
    x_linked = df['chrom'] == "chrX"
    y_linked = df['chrom'] == "chrY"

    early_onsets = df['HPO_IDs'].map(hpo_onset_modes)
    
    # For autosomal dominant disease, we can assign BS2 if the variant is observed in a healthy adult (either homozygous or heterozygous)
    autosomal_dominant = autosomal & dominant & ((df['gnomAD_joint_AF'] * df['gnomAD_joint_AN']) > 5) & np.logical_not(recessive) & np.logical_not(non_monogenic) & np.logical_not(non_mendelian)
    autosomal_recessive = autosomal & recessive & (((df['gnomAD_nhomalt_XX'] + df['gnomAD_nhomalt_XY']) > 5)) & np.logical_not(non_monogenic) & np.logical_not(non_mendelian)

    # For X-linked disease, we can assign BS2 if the variant is observed in a healthy adult male (hemizygous) or a healthy adult female (homozygous)
    x_linked_recessive = x_linked & recessive & (((df['gnomAD_nhomalt_XX'] + df['gnomAD_nhomalt_XY']) > 5)) & np.logical_not(non_monogenic) & np.logical_not(non_mendelian)
    x_linked_dominant = x_linked & dominant & np.logical_not(recessive) & (((df['gnomAD_nhomalt_XX'] + df['gnomAD_nhomalt_XY']) > 5) | ((df['gnomAD_joint_AF'] * df['gnomAD_joint_AN']) > 5)) & np.logical_not(non_monogenic) & np.logical_not(non_mendelian)

    y_linked = y_linked & (df['gnomAD_joint_AF_XY'] > 0) & np.logical_not(non_monogenic) & np.logical_not(non_mendelian)
    gnomad_bs2_criteria = autosomal_dominant | autosomal_recessive | x_linked_recessive | x_linked_dominant | y_linked
    bs2_criteria = gnomad_bs2_criteria & np.logical_not(bs1_criteria) & np.logical_not(incomplete_penetrance) & early_onsets
    return bs2_criteria




def BS4_criteria(df: pd.DataFrame, ped_df: pd.DataFrame, fam_name: str) -> pd.Series:
    # BS4: lack of family aggregation
    proband_info, father_info, mother_info, sib_info = identify_fam_members(ped_df, fam_name)
    proband, proband_pheno = proband_info
    father, father_pheno = father_info
    mother, mother_pheno = mother_info
    
    healthy_fam_members = []
    if father_pheno == 1:
        healthy_fam_members.append(father)
    if mother_pheno == 1:
        healthy_fam_members.append(mother)
    for sib, sib_pheno in sib_info.items():
        if sib_pheno == 1:
            healthy_fam_members.append(sib)

    final_criteria = np.array([False] * len(df))
    for healthy_mem in healthy_fam_members:
        final_criteria = final_criteria | (df[proband].str.count("1") <= df[healthy_mem].str.count("1"))

    return final_criteria


# Define this function at module level (outside any other function)
def process_gene_variants(args):
    """Helper function to unpack arguments for check_gene_variants"""
    return check_gene_variants(*args)


def check_gene_variants(gene, gene_df, pathogenic, proband, original_indices):
    """
    Check variants within a gene for trans/cis relationships with pathogenic variants.
    
    Args:
        gene: Gene name
        gene_df: DataFrame containing just variants for this gene
        pathogenic: Boolean array indicating which variants are pathogenic
        proband: Name of proband column
        original_indices: Array of indices in the original dataframe
        
    Returns:
        Tuple of arrays containing original indices where trans/cis conditions are true
    """
    pathogenic_variants = gene_df.loc[pathogenic, proband].tolist()
    logger.info(f"For gene {gene}, there are {len(pathogenic_variants)} pathogenic variants and {len(gene_df)} variants in total")
    
    var_in_trans = np.array([False] * len(gene_df))
    var_in_cis = np.array([False] * len(gene_df))
    
    if len([v for v in pathogenic_variants if len(v.split("|")) == 2 and v.split("|")[0] == "0"]) > 0:
        # If there is at least one pathogenic variant at the second copy of the proband's genome
        # Convert pandas Series to numpy array with np.array()
        var_in_trans = np.array((gene_df.loc[:, proband].str.split("|").str.get(0) == "1") & (gene_df.loc[:, "Gene"] == gene))
        var_in_cis = np.array((gene_df.loc[:, proband].str.split("|").str.get(1) == "1") & (gene_df.loc[:, "Gene"] == gene))
        logger.info(f"For gene {gene}, there are {var_in_trans.sum()} variants in-trans with pathogenic variants at the second copy of the proband's genome")
        logger.info(f"For gene {gene}, there are {var_in_cis.sum()} variants in-cis with pathogenic variants at the second copy of the proband's genome")

    if len([v for v in pathogenic_variants if len(v.split("|")) == 2 and v.split("|")[0] == "1"]) > 0:
        # If there is at least one pathogenic variant at the first copy of the proband's genome
        # Convert pandas Series to numpy array with np.array()
        var_in_trans_1 = np.array((gene_df.loc[:, proband].str.split("|").str.get(1) == "1") & (gene_df.loc[:, "Gene"] == gene))
        var_in_cis_1 = np.array((gene_df.loc[:, proband].str.split("|").str.get(0) == "1") & (gene_df.loc[:, "Gene"] == gene))
        logger.info(f"For gene {gene}, there are {var_in_trans_1.sum()} variants in-trans with pathogenic variants at the first copy of the proband's genome")
        logger.info(f"For gene {gene}, there are {var_in_cis_1.sum()} variants in-cis with pathogenic variants at the first copy of the proband's genome")
        var_in_trans = np.logical_or(var_in_trans, var_in_trans_1)
        var_in_cis = np.logical_or(var_in_cis, var_in_cis_1)

    # Map local boolean arrays to original indices
    trans_original_indices = original_indices[var_in_trans] if var_in_trans.any() else np.array([], dtype=int)
    cis_original_indices = original_indices[var_in_cis] if var_in_cis.any() else np.array([], dtype=int)
    logger.info(f"For gene {gene}, there are {len(trans_original_indices)} variants in-trans with pathogenic variants")
    logger.info(f"For gene {gene}, there are {len(cis_original_indices)} variants in-cis with pathogenic variants")
    
    return trans_original_indices, cis_original_indices


# Create a generator function that preserves gene order
def gene_args_generator(df, gene_to_indices, pathogenic, proband):
    for gene in df['Gene'].unique():
        indices = gene_to_indices[gene]
        yield (gene, 
                df.loc[indices, ["Gene", proband]], 
                pathogenic[indices], 
                proband, 
                indices)


def BP2_PM3_criteria(df: pd.DataFrame, 
                     ped_df: pd.DataFrame, 
                     fam_name: str,
                     gene_to_am_score_map: dict,
                     ps1_criteria: pd.Series,
                     ps2_criteria: pd.Series,
                     ps3_criteria: pd.Series,
                     is_recessive: np.ndarray,
                     is_dominant: np.ndarray,
                     is_non_monogenic: np.ndarray,
                     is_non_mendelian: np.ndarray,
                     haplo_insufficient: np.ndarray,
                     threads: int = 10) -> Tuple[pd.Series, pd.Series]:
    # BP2: observed in trans with a pathogenic variant in dominant disease, Or in-cis with a pathogenic variant with any inheritance mode
    # PM3: observed in trans with a pathogenic variant in recessive disease.
    pathogenic = df['vep_consq_lof'] | df['splicing_lof'] | df['5UTR_lof'] | ps1_criteria | (ps2_criteria & (df["CADD_phred"] >= 20)) | ps3_criteria

    # Groupby gene and see if some variant is in-trans or in-cis with a pathogenic variant using the pathogenic boolean array above
    proband_info, father_info, mother_info, sib_info = identify_fam_members(ped_df, fam_name)
    proband, proband_pheno = proband_info
    
    logger.info(f"Checking genes for in-trans and in-cis pathogenic variants")

    # Ensure Gene column is categorical for faster groupby
    df.loc[:, "Gene"] = df.loc[:, "Gene"].astype('category')
    
    # Create empty Series to hold results
    in_trans_pathogenic = pd.Series(False, index=df.index)
    in_cis_pathogenic = pd.Series(False, index=df.index)
    
    # Add pathogenic status to the dataframe temporarily
    df_with_pathogenic = df.assign(_pathogenic=pathogenic)
    
    # Define a function to process each gene group
    def process_gene_group(group):
        gene = group.name
        pathogenic_variants = group.loc[group['_pathogenic'], proband].tolist()
        logger.info(f"For gene {gene}, there are {len(pathogenic_variants)} pathogenic variants and {len(group)} variants in total")
        
        # Create Series to hold results for this gene
        trans_result = pd.Series(False, index=group.index)
        cis_result = pd.Series(False, index=group.index)
        
        # Check for variants on the second copy (where pathogenic has 0 in first position)
        path_second_copy = [v for v in pathogenic_variants if len(v.split("|")) == 2 and v.split("|")[0] == "0"]
        if len(path_second_copy) > 0:
            # Find variants in-trans (opposite haplotype) and in-cis (same haplotype)
            trans_mask = group[proband].str.split("|").str.get(0) == "1"
            cis_mask = group[proband].str.split("|").str.get(1) == "1"
            
            trans_result.loc[trans_mask] = True
            cis_result.loc[cis_mask] = True
            
            logger.info(f"For gene {gene}, there are {trans_mask.sum()} variants in-trans with pathogenic variants at the second copy")
            logger.info(f"For gene {gene}, there are {cis_mask.sum()} variants in-cis with pathogenic variants at the second copy")
        
        # Check for variants on the first copy (where pathogenic has 1 in first position)
        path_first_copy = [v for v in pathogenic_variants if len(v.split("|")) == 2 and v.split("|")[0] == "1"]
        if len(path_first_copy) > 0:
            # Find variants in-trans (opposite haplotype) and in-cis (same haplotype)
            trans_mask = group[proband].str.split("|").str.get(1) == "1"
            cis_mask = group[proband].str.split("|").str.get(0) == "1"
            
            # Update existing results with logical OR
            trans_result = trans_result | trans_mask
            cis_result = cis_result | cis_mask
            
            logger.info(f"For gene {gene}, there are {trans_mask.sum()} variants in-trans with pathogenic variants at the first copy")
            logger.info(f"For gene {gene}, there are {cis_mask.sum()} variants in-cis with pathogenic variants at the first copy")
        
        logger.info(f"For gene {gene}, there are {trans_result.sum()} variants in-trans with pathogenic variants")
        logger.info(f"For gene {gene}, there are {cis_result.sum()} variants in-cis with pathogenic variants")
        
        # Return a dataframe with two boolean columns for this gene group
        return pd.DataFrame({'trans': trans_result, 'cis': cis_result})
    
    # Apply the function to each gene group
    logger.info(f"Processing {df['Gene'].nunique()} genes with groupby.apply")
    results = df_with_pathogenic.groupby('Gene').apply(process_gene_group)
    
    # Unstack the results (convert MultiIndex to columns)
    if not results.empty:
        # Handle the case where results might have a MultiIndex
        if isinstance(results.index, pd.MultiIndex):
            # Combine the gene-level results back to the full dataframe
            for idx in results.index:
                gene_idx, variant_idx = idx
                in_trans_pathogenic.loc[variant_idx] = results.loc[idx, 'trans']
                in_cis_pathogenic.loc[variant_idx] = results.loc[idx, 'cis']
        else:
            # If only one gene was processed
            gene = results.index[0]
            in_trans_pathogenic.loc[results.index] = results['trans']
            in_cis_pathogenic.loc[results.index] = results['cis']
    
    # Convert Series back to numpy arrays before returning
    in_trans_pathogenic = in_trans_pathogenic.values
    in_cis_pathogenic = in_cis_pathogenic.values
    
    logger.info(f"There are {in_trans_pathogenic.sum()} variants in-trans with pathogenic variants")
    logger.info(f"There are {in_cis_pathogenic.sum()} variants in-cis with pathogenic variants")
    
    bp2_criteria = (in_trans_pathogenic & is_dominant) | in_cis_pathogenic
    pm3_criteria = in_trans_pathogenic & (is_recessive | is_non_monogenic)

    return bp2_criteria, pm3_criteria


def find_overlaps_bedtools_efficient(variants_df, regions_file):
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
    intersect_result = variants_bed.intersect(regions_bed, wa=True, u=True)
    
    # Extract the variant IDs that had overlaps
    overlapping_variants = set()
    for feature in intersect_result:
        # The variant_id is the 4th field (index 3)
        overlapping_variants.add(feature[3])
    
    return overlapping_variants


def BP3_criteria(df: pd.DataFrame, repeat_region_file: str, interpro_entry_map_pkl: str) -> pd.Series:
    # BP3: in-frame deletion in a repetitive region without a known function
    inframe_del = (df['Consequence'].str.contains('inframe_deletion')) | \
                  (df['Consequence'].str.contains('inframe_insertion')) | \
                  (df['vep_consq_length_changing'] & ~df['vep_consq_lof']) | \
                  (df['splicing_len_changing'] & ~df['splicing_lof']) | \
                  (df['5UTR_len_changing'] & ~df['5UTR_lof'])

    # Repeat region file is a gzipped bed file, we can read it with pandas
    df["variant_id"] = df["chrom"] + ":" + df["pos"].astype(str) + ":" + df["ref"] + "-" + df["alt"]
    in_repeat_regions = find_overlaps_bedtools_efficient(df, repeat_region_file)
    dm_instance = DomainNormalizer()
    interpro_entry_map_dict = pickle.load(open(interpro_entry_map_pkl, "rb"))
    functiona_domain = df.apply(lambda row:any(dm_instance.interpret_functionality(domain, interpro_entry_map_dict) == "Functional" for domain in str(row["DOMAINS"]).split("&")), axis=1)
    repetitive_region = ( df['DOMAINS'].str.contains('Low_complexity') | df['variant_id'].isin(in_repeat_regions) ) & ~functiona_domain
    not_deleterious = df['CADD_phred'] < 15
    return inframe_del & repetitive_region & not_deleterious


def BP5_criteria(df: pd.DataFrame, 
                 alt_disease_vcf: str, 
                 gene_to_am_score_map: dict, 
                 threads: int = 10) -> np.ndarray:
    '''
    BP5: variant found in a sample with known alternative molecular basis for disease
    Considers inheritance mode when determining presence in alternative disease samples
    '''
    # Convert DataFrame to dictionaries for parallel processing
    row_dicts = df.to_dict('records')
    
    # Prepare arguments for starmap
    args = [(row_dict, alt_disease_vcf, gene_to_am_score_map.get(row_dict['Gene'], np.nan)) for row_dict in row_dicts]
    
    # Process in parallel
    with mp.Pool(threads) as pool:
        results = pool.starmap(identify_presence_in_alt_disease_vcf, args)
    
    return np.array(results)
    

def identify_presence_in_alt_disease_vcf(row: dict, alt_disease_vcf: str, gene_mean_am_score: float) -> bool:
    '''
    Identify if the variant is present in the alt_disease_vcf with pysam.
    The alt_disease_vcf is a VCF file with known alternative molecular basis for disease.
    We need to identify whether the current variant is in a gene leading to recessive/dominant disease.
    
    For dominant disease:
        - Check if variant is present in alt_disease_vcf (heterozygous or homozygous)
        - If present, return True
    
    For recessive disease:
        - Check if variant is present with homozygous genotype in alt_disease_vcf
        - If present as homozygous, return True
    
    Args:
        row: Dictionary containing variant information
        alt_disease_vcf: Path to VCF file with known alternative molecular basis for disease
        gene_to_am_score_map: Dictionary mapping genes to AM scores for inheritance determination
        
    Returns:
        bool: True if variant is present according to inheritance mode criteria
    '''
    try:
        # First determine inheritance mode
        is_recessive, is_dominant, is_non_monogenic, is_non_mendelian, haplo_insufficient, incomplete_penetrance = identify_inheritance_mode_per_row(row, gene_mean_am_score)
        
        if not (is_dominant or is_recessive):
            is_dominant = False
            is_recessive = True
            
        # Open VCF file
        with pysam.VariantFile(alt_disease_vcf) as vcf:
            # Extract variant information from row
            chrom = row['chrom']  # Remove 'chr' prefix if present
            pos = int(row['pos']) # 0-based
            ref = row['ref']
            alt = row['alt']
            
            # Fetch variants in the region, remember that the region is 0-based and half-open, meaning the end is not included
            # End should be pos + length of ref (half-open leads to the exclusion of the end position)
            region = f"{chrom}:{max(1, pos-1)}-{pos+len(ref)}"
            
            for record in vcf.fetch(region=region):
                # Check exact position and allele match
                if record.pos == pos and record.ref == ref and alt in record.alts:
                    # For dominant inheritance, any presence is sufficient
                    if is_dominant:
                        if any(s.get('GT', (0,0)) != (0,0) for s in record.samples.values()):
                            return True
                    
                    # For recessive inheritance, need homozygous genotype
                    if is_recessive:
                        for sample in record.samples.values():
                            gt = sample.get('GT', (0,0))
                            # Check for homozygous alternate (1/1, 2/2, etc.)
                            if gt and len(gt) == 2 and gt[0] == gt[1] and gt[0] != 0:
                                return True
            
            return False
            
    except (ValueError, KeyError, OSError) as e:
        logger.warning(f"Error checking variant in alt disease VCF: {str(e)}")
        return False


def BP7_criteria(df: pd.DataFrame) -> pd.Series:
    # BP7: synonymous variant, no splicing-altering consequence, not conserved. 
    synonymous = df['Consequence'] == 'synonymous_variant'
    no_splicing_altering = (df['splicing_len_changing'] == False) | (df['splicing_lof'].isna()) | (df['splicing_lof'] == False) | (df['splicing_len_changing'].isna())
    not_conserved = df['Conservation'] <= 5 # 5 is the cutoff for highly conserved of a GERP++ score. 
    return synonymous & no_splicing_altering & not_conserved


def identify_fam_members(ped_df: pd.DataFrame, fam_name: str) -> pd.DataFrame:
    fam_ped_df = ped_df[ped_df['#FamilyID'] == fam_name]
    fam_ped_df["Phenotype"] = fam_ped_df["Phenotype"].astype(int)
    # fam_members = fam_ped_df['IndividualID'].tolist()
    father = fam_ped_df.loc[fam_ped_df["PaternalID"] != "0", "IndividualID"].tolist()
    mother = fam_ped_df.loc[fam_ped_df["MaternalID"] != "0", "IndividualID"].tolist()
    proband = fam_ped_df.loc[fam_ped_df["Phenotype"] == 2, "IndividualID"].tolist()[0]

    father = father[0] if father else None
    father_pheno = fam_ped_df.loc[fam_ped_df["IndividualID"] == father, "Phenotype"].tolist()[0] if father else None
    mother = mother[0] if mother else None
    mother_pheno = fam_ped_df.loc[fam_ped_df["IndividualID"] == mother, "Phenotype"].tolist()[0] if mother else None

    # Process siblings
    exist_mems = [m for m in [father, mother, proband] if m is not None]
    sibs = fam_ped_df.loc[~fam_ped_df["IndividualID"].isin(exist_mems), "IndividualID"].tolist()
    sib_pheno = fam_ped_df.loc[fam_ped_df["IndividualID"].isin(sibs), "Phenotype"].tolist()
    sib_info = dict(zip(sibs, sib_pheno))

    return (proband, 2), (father, father_pheno), (mother, mother_pheno), sib_info


def create_criteria_summary(row, criteria_order):
    """
    Create a summary string of active criteria for a given row.
    
    Args:
        row: A row from the criteria matrix
        criteria_order: List of criteria in order
        
    Returns:
        A semicolon-separated string of active criteria
    """
    active_criteria = [col for col in criteria_order if row[col]]
    return ";".join(active_criteria) if active_criteria else "None"



def summarize_acmg_criteria(df: pd.DataFrame, criteria_dict: Dict[str, np.ndarray]) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """
    Create a summary column of assigned criteria and a boolean matrix of all criteria.
    
    Args:
        df: Input DataFrame
        criteria_dict: Dictionary mapping criteria names to boolean arrays
        
    Returns:
        Tuple of:
        - Original DataFrame with added criteria summary column
        - Boolean DataFrame of all criteria (rows=variants, columns=criteria)
    """
    # Define criteria order
    criteria_order = [
        # Very Strong
        'PVS1',
        # Strong
        'PS1', 'PS2', 'PS3', 'PS4',
        # Moderate
        'PM1', 'PM2', 'PM3', 'PM4', 'PM5', 'PM6',
        # Supporting
        'PP1', 'PP2', 'PP3', 'PP4', 'PP5',
        # Stand-alone
        'BA1',
        # Strong Benign
        'BS1', 'BS2', 'BS3', 'BS4',
        # Supporting Benign
        'BP1', 'BP2', 'BP3', 'BP4', 'BP5', 'BP6', 'BP7'
    ]
    
    # Create criteria matrix
    criteria_matrix = pd.DataFrame(
        {name: criteria_dict.get(name, np.zeros(len(df), dtype=bool)) 
         for name in criteria_order},
        index=df.index
    )
    
    # Add summary column to original DataFrame
    df['ACMG_criteria'] = criteria_matrix.apply(create_criteria_summary, axis=1, criteria_order=criteria_order)
    
    return df, criteria_matrix



def calculate_posterior_probability(row, prior_probability=0.1, exp_base=2, odds_pvst=350):
    """
    Calculates the posterior probability of pathogenicity for a variant based on its criteria assignment using a Bayesian framework.

    Args:
        criteria_assignment (dict): A dictionary representing the strength of each evidence criterion.
                                     Keys are criteria names (e.g., "PVS1", "PS1", "PM2", "PP3", "BS1", "BP4"),
                                     and values are the corresponding strength scores.
        prior_probability (float): The prior probability of pathogenicity (default: 0.1).
        evidence_weights (dict): A dictionary mapping criteria names to their numerical weights, 
                                 according to the latest ClinGen SVI recommendations (default: illustrative values).
                                 PM2 is downgraded to a PP criterion

    Returns:
        float: The posterior probability of pathogenicity.
    """
    if isinstance(row, pd.Series):
        # Convert the Series to a dictionary
        criteria_assignment = row.to_dict()

    evidence_weights = {"PVS1": 1,
                        "PS1": 1/exp_base, "PS2": 1/exp_base, "PS3": 1/exp_base, "PS4": 1/exp_base,
                        "PM1": 1/exp_base**2, "PM2": 1/exp_base**3, "PM3": 1/exp_base**2, "PM4": 1/exp_base**2, "PM5": 1/exp_base**2, "PM6": 1/exp_base**2,
                        "PP1": 1/exp_base**3, "PP2": 1/exp_base**3, "PP3": 1/exp_base**3, "PP4": 1/exp_base**3, "PP5": 1/exp_base**3,
                        "BA1": -1,
                        "BS1": -1/exp_base, "BS2": -1/exp_base, "BS3": -1/exp_base, "BS4": -1/exp_base,
                        "BP1": -1/exp_base**3, "BP2": -1/exp_base**3, "BP3": -1/exp_base**3, "BP4": -1/exp_base**3, "BP5": -1/exp_base**3, "BP6": -1/exp_base**3, "BP7": -1/exp_base**3}

    # Default odds function (Illustrative - replace with a calibrated function based on SVI guidelines)
    odds_function = lambda total_weight: np.power(odds_pvst, total_weight)

    # There are 2 internal inconsistencies in the SVI guidelines.
    # Pathogenic(ii) and Likely Pathogenic(i) do not generate a posterior probability of the corresponding class.
    # We need to handle this by setting the odds of pathogenic to 0 if the total weight is negative. If we run into 2 PS evidences, we adjust the total_weight_pathogenic to 1514
    # If there is a PVS + 1 PM, we need to adjust their sum to 350

    ps_matches = [c for c,o in criteria_assignment.items() if c.startswith("PS") and int(o) > 0]
    ba_matches = [c for c,o in criteria_assignment.items() if c.startswith("BA") and int(o) > 0]
    pp_matches = [c for c,o in criteria_assignment.items() if (c.startswith("PP") or c.startswith("PM2")) and int(o) > 0]
    pm_matches = [c for c,o in criteria_assignment.items() if c.startswith("PM") and int(o) > 0 and not c.startswith("PM2")]
    pvs_matches = [c for c,o in criteria_assignment.items() if c.startswith("PVS") and int(o) > 0]
    bs_matches = [c for c,o in criteria_assignment.items() if c.startswith("BS") and int(o) > 0]
    bp_matches = [c for c,o in criteria_assignment.items() if c.startswith("BP") and int(o) > 0]

    odds_benign = odds_function(sum([evidence_weights[c] for c in bs_matches + bp_matches]))
    if len(ps_matches) >= 2:
        ps_odds = 1514 * odds_function(sum([evidence_weights[c] for c in ps_matches[:-2]]))
    else:
        ps_odds = odds_function(sum([evidence_weights[c] for c in ps_matches]))
    
    if len(pvs_matches) > 0 and len(pm_matches) > 0:
        pv_pm_odds = 350 * odds_function(sum([evidence_weights[c] for c in pm_matches[:-1]]))
    else:
        pv_pm_odds = odds_function(sum([evidence_weights[c] for c in pvs_matches + pm_matches]))

    pp_odds = odds_function(sum([evidence_weights[c] for c in pp_matches]))

    odds_path = ps_odds * pp_odds * pv_pm_odds

    # Calculate odds of pathogenicity
    odds_path = odds_path * odds_benign

    # Calculate posterior probability
    posterior_probability = (odds_path * prior_probability) / ((odds_path - 1) * prior_probability + 1)
    # Round posterior probability to 3 decimal places
    posterior_probability = round(posterior_probability, 3)

    if len(ba_matches) > 0:
        posterior_probability = 0

    # Classify the variant based on the posterior probability
    if posterior_probability >= 0.99:
        acmg_class = "Pathogenic"
    elif posterior_probability >= 0.9:
        acmg_class = "Likely Pathogenic"
    elif posterior_probability <= 0.001:
        acmg_class = "Benign"
    elif posterior_probability <= 0.1:
        acmg_class = "Likely Benign"
    else:
        acmg_class = "Uncertain Significance"

    return posterior_probability, acmg_class



# def calculate_posterior_probability(row, prior_probability=0.1, exp_base=2, odds_pvst=350):
#     """
#     Calculates the posterior probability of pathogenicity for a variant based on its criteria assignment using a Bayesian framework.

#     Args:
#         criteria_assignment (dict): A dictionary representing the strength of each evidence criterion.
#                                      Keys are criteria names (e.g., "PVS1", "PS1", "PM2", "PP3", "BS1", "BP4"),
#                                      and values are the corresponding strength scores.
#         prior_probability (float): The prior probability of pathogenicity (default: 0.1).
#         exp_base (float): The base of the exponential function (default: 2).
#         odds_pvst (float): The odds of pathogenicity (default: 350).

#     Returns:
#         float: The posterior probability of pathogenicity.
#     """
#     if isinstance(row, pd.Series):
#         # Convert the Series to a dictionary
#         criteria_assignment = row.to_dict()

#     evidence_weights = {"PVS1": 1,
#                         "PS1": 1/exp_base, "PS2": 1/exp_base, "PS3": 1/exp_base, "PS4": 1/exp_base,
#                         "PM1": 1/exp_base**2, "PM2": 1/exp_base**3, "PM3": 1/exp_base**2, "PM4": 1/exp_base**2, "PM5": 1/exp_base**2, "PM6": 1/exp_base**2,
#                         "PP1": 1/exp_base**3, "PP2": 1/exp_base**3, "PP3": 1/exp_base**3, "PP4": 1/exp_base**3, "PP5": 1/exp_base**3,
#                         "BA1": -1,
#                         "BS1": -1/exp_base, "BS2": -1/exp_base, "BS3": -1/exp_base, "BS4": -1/exp_base,
#                         "BP1": -1/exp_base**3, "BP2": -1/exp_base**3, "BP3": -1/exp_base**3, "BP4": -1/exp_base**3, "BP5": -1/exp_base**3, "BP6": -1/exp_base**3, "BP7": -1/exp_base**3}

#     # Default odds function (Illustrative - replace with a calibrated function based on SVI guidelines)
#     odds_function = lambda total_weight: np.power(odds_pvst, total_weight)

#     # Calculate total weight based on provided evidence weights
#     total_weight_pathogenic = 0
#     total_weight_benign = 0
#     for criterion, occurence in criteria_assignment.items():
#         if criterion == "BA1" and occurence > 0:
#             total_weight_pathogenic = 0
#             total_weight_benign = 0
#             break
#         if criterion in evidence_weights:
#             occurence = 1 if int(occurence) > 0 else 0
#             if criterion.startswith("P"):  # Pathogenic criteria
#                 total_weight_pathogenic += evidence_weights[criterion] * occurence
#             elif criterion.startswith("B"):  # Benign criteria
#                 total_weight_benign += evidence_weights[criterion] * occurence
#         else:
#             continue

#     # Calculate odds of pathogenicity
#     odds_path = odds_function(total_weight_pathogenic + total_weight_benign)

#     # Calculate posterior probability
#     posterior_probability = (odds_path * prior_probability) / ((odds_path - 1) * prior_probability + 1)

#     # Classify the variant based on the posterior probability
#     if posterior_probability >= 0.99:
#         acmg_class = "Pathogenic"
#     elif posterior_probability >= 0.90:
#         acmg_class = "Likely Pathogenic"
#     elif posterior_probability <= 0.001:
#         acmg_class = "Benign"
#     elif posterior_probability <= 0.1:
#         acmg_class = "Likely Benign"
#     else:
#         acmg_class = "Uncertain Significance"

#     return posterior_probability, acmg_class



def sort_and_rank_variants(df: pd.DataFrame) -> pd.DataFrame:
    """
    Sort variants by their maximum ACMG quantitative score and add ranking.
    
    Args:
        df: DataFrame with ACMG_quant_score column and variant information
        
    Returns:
        DataFrame sorted by variant's max ACMG score with added variant rank column
    """
    # Group by variant coordinates to get max score per variant
    variant_groups = df.groupby(['chrom', 'pos', 'ref', 'alt'])
    max_scores = variant_groups['ACMG_quant_score'].transform('max')
    
    # Sort the DataFrame by max scores (descending) while maintaining variant grouping
    df['max_variant_score'] = max_scores
    df_sorted = df.sort_values(
        by=['max_variant_score', 'chrom', 'pos', 'ref', 'alt', 'ACMG_quant_score'],
        ascending=[False, True, True, True, True, False]
    )
    
    # Add variant rank (same rank for all rows of same variant)
    variant_rank = (df_sorted.groupby(['chrom', 'pos', 'ref', 'alt'])
                            .ngroup()
                            .reset_index()
                            .set_index('index') + 1)
    df_sorted['variant_rank'] = variant_rank
    
    # Clean up temporary column
    df_sorted = df_sorted.drop('max_variant_score', axis=1)
    
    return df_sorted


def BA1_criteria(anno_df: pd.DataFrame) -> pd.Series:
    """
    Apply BA1 criteria to the annotation DataFrame.
    """
    false_neg_rate, common_vars = control_false_neg_rate(anno_df['gnomAD_joint_AF_max'], anno_df['gnomAD_joint_AN_max'], af_threshold=0.05, alpha=0.01)
    ba1_criteria = np.where(common_vars.isna(), anno_df['gnomAD_joint_AF'] > 0.05, common_vars)
    logger.info(f"BA1 criteria applied, {ba1_criteria.sum()} variants are having the BA1 criteria")

    return ba1_criteria



def ACMG_criteria_assign(anno_table: str, 
                         am_score_table: str, 
                         clinvar_patho_af_stat: str,
                         clinvar_patho_exon_af_stat: str,
                         clinvar_aa_dict_pkl: str,
                         clinvar_splice_dict_pkl: str,
                         interpro_entry_map_pkl: str,
                         intolerant_domains_pkl: str,
                         intolerant_motifs_pkl: str,
                         domain_mechanism_tsv: str,
                         tranx_exon_domain_map_pkl: str,
                         repeat_region_file: str,
                         mavedb_metadata_tsv: str = "",
                         fam_name: str = "",
                         ped_table: str = "",
                         alt_disease_vcf: str = "",
                         gnomAD_extreme_rare_threshold: float = 0.0001,
                         expected_incidence: float = 0.001,
                         threads: int = 10) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """
    Main function to assign ACMG criteria.
    Returns both annotated DataFrame and criteria matrix.

    The ClinGen SVI group has introduced a new dimension of evaluation for each criteria called Strength of Criteria.
    But many strengths cannot be applied at bioinformatic level. 

    ===========Regarding the combining rule===========
    Regarding the combining rule, SVI from ClinGen suggests that PVS1 + 1PP = Likely Pathogenic.
    We just adopted the naive Bayesian Framework to calculate the posterior probability of pathogenicity.

    ===========Regarding the criteria, Cannot be applied===========
    1. PS2/PM6, regarding the DeNovo related evidence. The strength is related to the phenotype consistency with the disease. Such clinical information is rarely normalized for automatic interpretation. Thus Skip.
    2. PP1, the refinement not applicable because multiple family co-segregation information is rarely available for a same type of disease. 
    3. PP4, Patient's phenotype or family history is highly specific for a disease with a single genetic etiology. (This cannot be applied simutaneously with PP1, because if a disease is tightly linked with only one gene, the segregation is doomed. Therefore further segregation observation does not add any more confidence because the confidence from this perspective is already reaching a ceiling.)
    4. PS3/BS3, the functional assay for a specific variant is so rare and so hard to fetch, nearly impossible for practice.

    ===========Regarding the criteria, Can be applied===========
    1. PVS1, LOEUF <= 0.35 should be considered as intolerant to LoF.
    2. PM2, is reduced from Moderate to Supporting.
    3. PP5/BP6, reputable source reported as benign or pathogenic. Suggested to be abandoned or at least not assigned along with PS3/BS3.
    3. PM3, can be only applied if PM2 is True (sufficiently rare in gnomAD).
    """
    anno_df = pd.read_table(anno_table, low_memory=False)
    logger.info(f"Got {threads} threads to process the input table {anno_table}, now the table looks like: \n{anno_df[:5].to_string(index=False)}")
    anno_df = vep_consq_interpret(anno_df, threads)

    am_score_df = pd.read_table(am_score_table, low_memory=False)
    if ped_table:
        ped_df = pd.read_table(ped_table, low_memory=False)
    else:
        ped_df = None

    # Convert the am_score_df to a dictionary:
    # 1. Ensembl transcript ID (column 'transcript') to mean AM score (column 'mean_am_pathogenicity')
    am_score_dict = dict(zip(am_score_df['transcript'], am_score_df['mean_am_pathogenicity']))
    # Use anno_df to create a dict map from Ensembl transcript ID to gene ID
    transcript_to_gene_map = dict(zip(anno_df['Feature'], anno_df['Gene']))
    # Use the two dict above to create dict that maps gene ID to mean AM score
    gene_to_am_score_map = {g: am_score_dict[t] for t, g in transcript_to_gene_map.items() if t in am_score_dict}
    clinvar_aa_dict = pickle.load(open(clinvar_aa_dict_pkl, "rb"))
    clinvar_aa_gene_map = {g: clinvar_aa_dict[t] for t, g in transcript_to_gene_map.items() if t in clinvar_aa_dict}
    logger.info(f"gene_to_am_score_map created, {len(gene_to_am_score_map)} genes are having the AM score")

    # Apply the PVS1 criteria, LoF on a gene known to to be pathogenic due to LoF
    pvs1_criteria = PVS1_criteria(anno_df, 
                                  gene_to_am_score_map, 
                                  clinvar_aa_gene_map, 
                                  clinvar_patho_exon_af_stat,
                                  interpro_entry_map_pkl,
                                  ped_df,
                                  tranx_exon_domain_map_pkl=tranx_exon_domain_map_pkl) # When test on ClinVar variants, fam_name is set to None because no genotype information are provided
    logger.info(f"PVS1 criteria applied, {pvs1_criteria.sum()} variants are having the PVS1 criteria")
    gc.collect()

    # Apply the PS3 and BS3 criteria
    ps3bs3_results, anno_df = PS3_BS3_criteria(anno_df, pvs1_criteria, mavedb_metadata_tsv)
    logger.info(f"PS3 criteria applied, {ps3bs3_results['PS3'].sum()} variants are having the PS3 criteria")

    # Prevent double counting, if PVS1 is True, then PS3 should be False
    ps3_criteria = ps3bs3_results['PS3']
    bs3_criteria = ps3bs3_results['BS3']
    logger.info(f"PS3 criteria applied, after preventing double counting with PVS1, {ps3_criteria.sum()} variants are having the PS3 criteria")
    logger.info(f"BS3 criteria applied, {bs3_criteria.sum()} variants are having the BS3 criteria")
    gc.collect()

    # Apply the PS1 and PM5 criteria
    ps1_criteria, pm5_criteria = PS1_PM5_criteria(anno_df, clinvar_aa_dict_pkl, clinvar_splice_dict_pkl, threads)
    pm5_criteria = pm5_criteria & ~ps1_criteria
    ps1_criteria = ps1_criteria & ~pvs1_criteria & ~ps3_criteria
    logger.info(f"PS1 criteria applied, {ps1_criteria.sum()} variants are having the PS1 criteria") 
    logger.info(f"PM5 criteria applied, {pm5_criteria.sum()} variants are having the PM5 criteria")
    gc.collect()
    # Apply the PS2 and PM6 criteria
    if not ped_df is None and not fam_name is None:
        ps2_criteria, pm6_criteria = PS2_PM6_criteria(anno_df, ped_df, fam_name, threads=threads)
    else:
        logger.warning(f"No ped_table provided, skip the PS2 and PM6 criteria")
        ps2_criteria, pm6_criteria = np.array([False] * len(anno_df)), np.array([False] * len(anno_df))
    logger.info(f"PS2 criteria applied, {ps2_criteria.sum()} variants are having the PS2 criteria")
    logger.info(f"PM6 criteria applied, {pm6_criteria.sum()} variants are having the PM6 criteria")
    gc.collect()
   
    '''
    PS4 cannot be applied because usually we dont have enough cases to determine the frequency of the variant
    '''
    # Apply PM1 criteria, mutational hotspot or well-established functional protein domain
    pm1_criteria = PM1_criteria(anno_df, intolerant_domains_pkl, intolerant_motifs_pkl, threads)
    pm1_criteria = pm1_criteria & ~pvs1_criteria & ~ps1_criteria
    logger.info(f"PM1 criteria applied, {pm1_criteria.sum()} variants are having the PM1 criteria")
    gc.collect()
    # Apply PM2 criteria, absent from gnomAD or extremely rare in gnomAD
    pm2_criteria = PM2_criteria(anno_df, gnomAD_extreme_rare_threshold)
    logger.info(f"PM2 criteria applied, {pm2_criteria.sum()} variants are having the PM2 criteria")
    gc.collect()
    # Apply PM4 criteria, causing the protein length change
    pm4_criteria = PM4_criteria(anno_df)
    gc.collect()
    # Prevent double counting of PM4
    pm4_criteria = pm4_criteria & ~pvs1_criteria
    logger.info(f"PM4 criteria applied, {pm4_criteria.sum()} variants are having the PM4 criteria")
    '''
    PP1 cannot be applied because usually we dont have multiple family segregation information for each variant
    '''
    # Apply PP2 criteria, missense variant in a gene/domain that not only intolerant to truncating variants but also intolerant to missense variants
    pp2_criteria, bp1_criteria = PP2_BP1_criteria(anno_df, domain_mechanism_tsv)
    logger.info(f"PP2 criteria applied, {pp2_criteria.sum()} variants are having the PP2 criteria")
    logger.info(f"BP1 criteria applied, {bp1_criteria.sum()} variants are having the BP1 criteria")
    gc.collect()
    # Apply PP3 criteria, predicted to be deleterious by in-silico tools
    pp3_criteria, bp4_criteria = PP3_BP4_criteria(anno_df)
    bp4_criteria = bp4_criteria & ~bs3_criteria
    logger.info(f"BP4 criteria applied, {bp4_criteria.sum()} variants are having the BP4 criteria")
    gc.collect()
    # Prevent double counting of PP3
    pp3_criteria = pp3_criteria & ~pvs1_criteria & ~pm4_criteria
    logger.info(f"PP3 criteria applied, {pp3_criteria.sum()} variants are having the PP3 criteria")
    gc.collect()
    '''
    PP4 cannot be applied, Patient's phenotype or family history is highly specific for a disease with a single genetic etiology
    '''
    # Apply PP5 criteria, reported as pathogenic by a reputable source but without to many supporting evidences
    # Apply BP6 criteria, reported as benign by a reputable source but without to many supporting evidences
    pp5_criteria, bp6_criteria = PP5_BP6_criteria(anno_df)
    pp5_criteria = pp5_criteria & ~ps3_criteria
    bp6_criteria = bp6_criteria & ~bs3_criteria
    logger.info(f"PP5 criteria applied, {pp5_criteria.sum()} variants are having the PP5 criteria")
    logger.info(f"BP6 criteria applied, {bp6_criteria.sum()} variants are having the BP6 criteria")
    gc.collect()
    
    ba1_criteria = BA1_criteria(anno_df)

    # Apply BS1, PAF of variant is greater than expected incidence of the disease
    bs1_criteria, recessive, dominant, non_monogenic, non_mendelian, haplo_insufficient, incomplete_penetrance = BS1_criteria(anno_df, gene_to_am_score_map, threads = threads, expected_incidence = expected_incidence, clinvar_patho_af_stat = clinvar_patho_af_stat)
    logger.info(f"BS1 criteria applied, {bs1_criteria.sum()} variants are having the BS1 criteria")
    gc.collect()

    # Apply BS2, variant observed in a healthy adult
    bs2_criteria = BS2_criteria(anno_df, bs1_criteria, recessive, dominant, non_monogenic, non_mendelian, haplo_insufficient, incomplete_penetrance)
    logger.info(f"BS2 criteria applied, {bs2_criteria.sum()} variants are having the BS2 criteria")
    gc.collect()

    # Apply BS4, lack of family segregation
    if not ped_df is None and not fam_name is None:
        bs4_criteria = BS4_criteria(anno_df, ped_df, fam_name)
    else:
        logger.warning(f"No ped_table provided, skip the BS4 criteria")
        bs4_criteria = np.array([False] * len(anno_df))
    logger.info(f"BS4 criteria applied, {bs4_criteria.sum()} variants are having the BS4 criteria")
    # Apply BP3, in-frame indels in a repetitive region without a known function
    bp3_criteria = BP3_criteria(anno_df, repeat_region_file, interpro_entry_map_pkl)
    logger.info(f"BP3 criteria applied, {bp3_criteria.sum()} variants are having the BP3 criteria")
    gc.collect()
    # Apply BP5, variant found in a sample with known alternative molecular basis for disease
    if alt_disease_vcf:
        bp5_criteria = BP5_criteria(anno_df, alt_disease_vcf, gene_to_am_score_map, threads)
    else:
        bp5_criteria = np.array([False] * len(anno_df))
    logger.info(f"BP5 criteria applied, {bp5_criteria.sum()} variants are having the BP5 criteria")
    gc.collect()
    # Apply BP7, synonymous variant, no splicing-altering consequence, not conserved. 
    bp7_criteria = BP7_criteria(anno_df)
    logger.info(f"BP7 criteria applied, {bp7_criteria.sum()} variants are having the BP7 criteria")
    gc.collect()
    # Apply BP2, observed in trans with a pathogenic variant in dominant disease, Or in-cis with a pathogenic variant with any inheritance mode
    # Apply PM3, observed in trans with a pathogenic variant in recessive disease.
    if not ped_df is None and not fam_name is None:
        bp2_criteria, pm3_criteria = BP2_PM3_criteria(anno_df, 
                                                      ped_df, 
                                                      fam_name,
                                                      gene_to_am_score_map,
                                                      ps1_criteria,
                                                      ps2_criteria,
                                                      ps3_criteria,
                                                      recessive,
                                                      dominant,
                                                      non_monogenic,
                                                      non_mendelian,
                                                      haplo_insufficient,
                                                      threads)
    else:
        logger.warning(f"No ped_table provided, skip the BP2 and PM3 criteria")
        bp2_criteria, pm3_criteria = np.array([False] * len(anno_df)), np.array([False] * len(anno_df))
    logger.info(f"BP2 criteria applied, {bp2_criteria.sum()} variants are having the BP2 criteria")
    logger.info(f"PM3 criteria applied, {pm3_criteria.sum()} variants are having the PM3 criteria")
    gc.collect()
    # Collect all criteria in a dictionary
    criteria_dict = {
        'PVS1': pvs1_criteria,
        'PS1': ps1_criteria,
        'PS2': ps2_criteria,
        'PS3': ps3_criteria,
        'PM1': pm1_criteria,
        'PM2': pm2_criteria,
        'PM3': pm3_criteria,
        'PM4': pm4_criteria,
        'PM5': pm5_criteria,
        'PM6': pm6_criteria,
        'PP2': pp2_criteria,
        'PP3': pp3_criteria,
        'PP5': pp5_criteria,
        'BA1': ba1_criteria,
        'BS1': bs1_criteria,
        'BS2': bs2_criteria,
        'BS3': bs3_criteria,
        'BS4': bs4_criteria,
        'BP1': bp1_criteria,
        'BP2': bp2_criteria,
        'BP3': bp3_criteria,
        'BP4': bp4_criteria,
        'BP5': bp5_criteria,
        'BP6': bp6_criteria,
        'BP7': bp7_criteria
    }
    
    # Create summary and matrix
    anno_df, criteria_matrix = summarize_acmg_criteria(anno_df, criteria_dict)
    
    # Apply quantification using ACMG_criteria column and use that to sort the variants
    posterior_probability, acmg_class = zip(*criteria_matrix.apply(calculate_posterior_probability, axis=1))
    anno_df["ACMG_quant_score"] = posterior_probability
    anno_df["ACMG_class"] = acmg_class

    # Sort and rank variants
    anno_df = sort_and_rank_variants(anno_df)
    # Save the annotated table to replace the input anno_table
    anno_df.to_csv(anno_table, sep="\t", index=False)

    # Save the criteria matrix to a file which the path is based on the input anno_table
    output_matrix = ".".join(anno_table.split(".")[:-1]) + ".acmg.tsv"
    criteria_matrix.to_csv(output_matrix, sep="\t", index=False)
    
    return anno_df, criteria_matrix



if __name__ == "__main__":
    parser = ap.ArgumentParser()
    parser.add_argument("--anno_table", type=str, required=True)
    parser.add_argument("--am_score_table", type=str, required=True)
    parser.add_argument("--ped_table", type=str, required=False, default=None)
    parser.add_argument("--fam_name", type=str, required=False, default=None)
    parser.add_argument("--clinvar_patho_af_stat", type=str, required=True)
    parser.add_argument("--clinvar_patho_exon_af_stat", type=str, required=True)
    parser.add_argument("--clinvar_aa_dict_pkl", type=str, required=True)
    parser.add_argument("--clinvar_splice_dict_pkl", type=str, required=True)
    parser.add_argument("--interpro_entry_map_pkl", type=str, required=True)
    parser.add_argument("--intolerant_domains_pkl", type=str, required=True)
    parser.add_argument("--intolerant_motifs_pkl", type=str, required=True)
    parser.add_argument("--domain_mechanism_tsv", type=str, required=True)
    parser.add_argument("--tranx_exon_domain_map_pkl", type=str, required=True)
    parser.add_argument("--am_score_vcf", type=str, required=False, default=None)
    parser.add_argument("--alt_disease_vcf", type=str, required=False, default=None)
    parser.add_argument("--mavedb_metadata_tsv", type=str, required=False, default=None)
    parser.add_argument("--repeat_region_file", type=str, required=True)
    parser.add_argument("--gnomAD_extreme_rare_threshold", type=float, required=False, default=0.0001)
    parser.add_argument("--expected_incidence", type=float, required=False, default=0.001)
    parser.add_argument("--threads", type=int, required=False, default=10)
    args = parser.parse_args()

    anno_df, criteria_matrix = ACMG_criteria_assign(args.anno_table, 
                                                    args.am_score_table, 
                                                    args.clinvar_patho_af_stat,
                                                    args.clinvar_patho_exon_af_stat,
                                                    args.clinvar_aa_dict_pkl,
                                                    args.clinvar_splice_dict_pkl,
                                                    args.interpro_entry_map_pkl,
                                                    args.intolerant_domains_pkl,
                                                    args.intolerant_motifs_pkl,
                                                    args.domain_mechanism_tsv,
                                                    args.tranx_exon_domain_map_pkl,
                                                    args.repeat_region_file,
                                                    mavedb_metadata_tsv=args.mavedb_metadata_tsv,
                                                    fam_name=args.fam_name,
                                                    ped_table=args.ped_table,
                                                    alt_disease_vcf=args.alt_disease_vcf,
                                                    gnomAD_extreme_rare_threshold=args.gnomAD_extreme_rare_threshold,
                                                    expected_incidence=args.expected_incidence,
                                                    threads=args.threads)




