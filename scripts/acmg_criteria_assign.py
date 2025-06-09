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
from scipy import stats
from scipy.stats import binom
import mmap
import gc
import gzip

from stat_protein_domain_amscores import nested_defaultdict
from protein_domain_mapping import DomainNormalizer
from mavdb_interpreter import MaveDBScoreInterpreter
from find_cosegregation_vars import find_cosegregating_variants
from determine_phase import batch_annotate_cis_trans_from_table
from splicing_var_analysis import parse_hgvsc_splice_position

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
    nmd_escaping = "escaping" in str(row.get('NMD', ""))
    if not isinstance(consq, str):
        return False, False
        
    # Skip splicing related consequences to leave them handled by SpliceAI and SpliceVault
    lof_criteria = {
        'stop_gained', 
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
    
    is_lof = any(c in consq for c in lof_criteria) & (loftee_result != 'LC') & (nmd_escaping == False)

    # For pure splicing related consequences, we trust SpliceAI and SpliceVault more than VEP
    is_length_changing = is_lof or any(c in consq for c in length_changing_criteria) or any(c in consq for c in lof_criteria) or nmd_escaping
    
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



def summarize_clinvar_gene_pathogenicity(clinvar_gene_aa_dict: dict, high_confidence_status = {
        # Higher confidence (2+ stars)
        'practice_guideline': 4,                                   # 4 stars
        'reviewed_by_expert_panel': 3,                             # 3 stars
        'criteria_provided,_multiple_submitters,_no_conflicts': 2,  # 2 stars
    }) -> set:
    '''
    Based on the clinvar_gene_aa_dict:
    {ensg: {protein_pos: {hgvsp: {'CLNSIG': [cln_sig], 'CLNREVSTAT': [rev_stat]}}}}

    Extract the gene list that has been reported to harbor pathogenic variants in ClinVar database,
    considering both CLNSIG and CLNREVSTAT values.
    '''

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
        elif (np.logical_not(df["Consequence"].str.contains("start_lost")) & df["BIOTYPE"].fillna(".").str.contains("protein_coding")).any():
            return df["Gene"].values[0] if len(df["Gene"].values) > 0 else np.nan
        else:
            return np.nan
    else:
        return np.nan
    

def downstream_domain_impact(exon_str, tranx_id, tranx_exon_domain_map, interpro_entry_map_dict, dm_instance, domains="", intol_domains=[], ensg_id=""):
    '''
    Used for frameshift and stopgain variants to explore whether the downstream protein region involving functional domains
    '''
    affected_exons = set([])
    if not isinstance(exon_str, str):
        pass
    elif "-" in exon_str and "/" in exon_str:
        # e.g. 2-3/5
        affected_exons.update(range(min(int(exon_str.split("-")[1].split("/")[0])+1, int(exon_str.split("/")[1])), int(exon_str.split("/")[1])))
    elif "/" in exon_str:
        # e.g. 2/5
        affected_exons.update(range(min(int(exon_str.split("/")[0])+1, int(exon_str.split("/")[1])), int(exon_str.split("/")[1])))
    else:
        raise ValueError(f"Invalid exon string: {exon_str}")
    
    if isinstance(domains, str):
        domains = domains.split("&")
        for domain in domains:
            if ensg_id + ":" + domain in intol_domains:
                return True
            elif dm_instance.interpret_functionality(domain, interpro_entry_map_dict) == "Functional":
                return True
    
    tranx_id = tranx_id.split(".", 1)[0] # Remove the ENSG version number
    for exon in affected_exons:
        if tranx_id in tranx_exon_domain_map:
            if exon in tranx_exon_domain_map[tranx_id]:
                domains = tranx_exon_domain_map[tranx_id][exon]
                for domain in domains:
                    domain = domain.split(":", 1)[1] # Remove the ENSG prefix in the domain path
                    if ensg_id + ":" + domain in intol_domains:
                        return True
                    elif dm_instance.interpret_functionality(domain, interpro_entry_map_dict) == "Functional":
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
        affected_exons_patho_af.add(clinvar_patho_exon_af_dict.get(tranx_id, {}).get(exon, (0, ))[0])
    
    if logic == "any":
        return any(float(epa) < threshold for epa in affected_exons_patho_af)
    elif logic == "all":
        return all(float(epa) < threshold for epa in affected_exons_patho_af)
    else:
        raise ValueError(f"Invalid logic: {logic}, it should be either 'any' or 'all', depending on your needs")


def span_functional_domains(row, 
                            dm_instance=None, 
                            interpro_entry_map_dict=None, 
                            tranx_exon_domain_map=None, 
                            clinvar_patho_exon_af_dict=None, 
                            intol_domains=[],
                            exon_patho_af_threshold=0.01):
    '''
    Identify whether a truncating variant is involving functional regions on a protein
    '''
    ensg_id = str(row['Gene'])
    func_domain = False
    exon_frequent_patho = False
    if ("frameshift" in row['Consequence']) or ("stop_gained" in row['Consequence']):
        # These variants not only affect the local region of proteins, but also affect the downstream protein regions
        func_domain = downstream_domain_impact(row['EXON'], row['Feature'], tranx_exon_domain_map, interpro_entry_map_dict, dm_instance, domains=row["DOMAINS"], intol_domains=intol_domains, ensg_id=ensg_id)
        exon_frequent_patho = downstream_exon_patho_af(row, clinvar_patho_exon_af_dict, logic="any", threshold=exon_patho_af_threshold)
        logger.info(f"For variant {row['chrom']}:{row['pos']} with transcript {row['Feature']}, the exon is: {row['EXON']}, and the pathogenic AF is: {exon_frequent_patho}")
    elif row["5UTR_lof"] or row["5UTR_frameshift"] or row["5UTR_len_changing"]:
        func_domain = row["5UTR_span_intol_domain"]
        affected_exons = {"1"}
        affected_exon_patho_af = { exon: clinvar_patho_exon_af_dict.get(row['Feature'], {}).get(exon, (0, ))[0] for exon in affected_exons }
        logger.info(f"For variant {row['chrom']}:{row['pos']} with transcript {row['Feature']}, the affected exons are: {affected_exons}, and the pathogenic AFs are: {affected_exon_patho_af}")
        exon_frequent_patho = any(float(epa) < exon_patho_af_threshold for epa in affected_exon_patho_af.values())
    elif row["splicing_lof"] or row["splicing_frameshift"] or row["splicing_len_changing"]:
        func_domain = row["splicing_span_intol_domain"]
        affected_exons = set(str(row["splicing_affected_exons"]).split(","))
        affected_exon_patho_af = { exon: clinvar_patho_exon_af_dict.get(row['Feature'], {}).get(exon, (0, ))[0] for exon in affected_exons }
        logger.info(f"For variant {row['chrom']}:{row['pos']} with transcript {row['Feature']}, the affected exons are: {affected_exons}, and the pathogenic AFs are: {affected_exon_patho_af}")
        exon_frequent_patho = any(float(epa) < exon_patho_af_threshold for epa in affected_exon_patho_af.values())
    else:
        if isinstance(row["DOMAINS"], str):
            if intol_domains:
                func_domain = any([ensg_id + ":" + domain in intol_domains for domain in row["DOMAINS"].split("&")])
            else:
                func_domain = any([dm_instance.interpret_functionality(domain, interpro_entry_map_dict) == "Functional" for domain in row["DOMAINS"].split("&")])
        exon_frequent_patho = float(clinvar_patho_exon_af_dict.get(row['Feature'], {}).get(row['EXON'], (0, ))[0]) < exon_patho_af_threshold
        logger.debug(f"For variant {row['chrom']}:{row['pos']} with transcript {row['Feature']}, the exon is: {row['EXON']}, and the pathogenic AF is: {exon_frequent_patho}")
    
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
    total_protein_len = pd.to_numeric(df["Protein_position"].str.split("/", expand=True).iloc[:, 1], errors='coerce')
    truncate_frac = np.where(nonsense_vars, 1 - (var_protein_pos / total_protein_len), var_span / total_protein_len)
    return truncate_frac



def PVS1_criteria(df: pd.DataFrame, 
                  clinvar_gene_aa_dict: dict,
                  clinvar_patho_exon_af_stat: str,
                  interpro_entry_map_pkl: str,
                  gene_to_am_score_map: dict,
                  intolerant_domains: set = [],
                  tranx_exon_domain_map_pkl: str = None,
                  proband_gt_col: str = None) -> pd.DataFrame:
    '''
    Introducing varied strength levels of PVS1,
    PVS1: 4
    PVS1_strong: 3
    PVS1_moderate: 2
    PVS1_supporting: 1
    not_applicable: 0
    '''
    # Load the domains
    clinvar_pathogenic_genes = summarize_clinvar_gene_pathogenicity(clinvar_gene_aa_dict)
    clinvar_pathogenic = df['Gene'].isin(clinvar_pathogenic_genes)
    logger.info(f"{clinvar_pathogenic.sum()} variants are having ClinVar pathogenic variants")
    df["Gene_avg_AM_score"] = df['Gene'].map(gene_to_am_score_map)

    if proband_gt_col:
        heterozygous = df[proband_gt_col].str.count("1") == 1
        homozygous = df[proband_gt_col].str.count("1") == 2
    else:
        heterozygous = np.array([True] * len(df))
        homozygous = np.array([False] * len(df))
    
    lof_intol_metric = (df["LOEUF"].fillna(2) < 0.35) | (df["Gene_avg_AM_score"].fillna(0) > 0.7)
    logger.info(f"For LOEUF < 0.35, {(df["LOEUF"].fillna(2) < 0.35).sum()} variants are located in a gene intolerant to LoF variants")
    logger.info(f"For mean AM score > 0.7, {(df["Gene_avg_AM_score"].fillna(0) > 0.7).sum()} variants are located in a gene intolerant to LoF variants")
    lof_intolerant_het = (clinvar_pathogenic | lof_intol_metric) & heterozygous
    lof_intolerant_hom = clinvar_pathogenic & homozygous
    lof_intolerant = lof_intolerant_het | lof_intolerant_hom
    
    # Load the necessary dict file
    clinvar_patho_exon_af_dict = pickle.load(gzip.open(clinvar_patho_exon_af_stat)) if clinvar_patho_exon_af_stat.endswith(".gz") else pickle.load(open(clinvar_patho_exon_af_stat, 'rb'))
    interpro_entry_map_dict = pickle.load(gzip.open(interpro_entry_map_pkl)) if interpro_entry_map_pkl.endswith(".gz") else pickle.load(open(interpro_entry_map_pkl, 'rb'))
    tranx_exon_domain_map = pickle.load(gzip.open(tranx_exon_domain_map_pkl)) if tranx_exon_domain_map_pkl.endswith(".gz") else pickle.load(open(tranx_exon_domain_map_pkl, 'rb'))
    dm_instance = DomainNormalizer()
    intolerant_domains, exon_rare_patho_afs = zip(*df.apply(span_functional_domains, axis=1, dm_instance=dm_instance, 
                                                                                            interpro_entry_map_dict=interpro_entry_map_dict, 
                                                                                            tranx_exon_domain_map=tranx_exon_domain_map, 
                                                                                            clinvar_patho_exon_af_dict=clinvar_patho_exon_af_dict, 
                                                                                            intol_domains=intolerant_domains,
                                                                                            exon_patho_af_threshold=0.01))
    
    # Convert tuples to numpy arrays for boolean operations
    intolerant_domains = np.array(intolerant_domains)
    exon_rare_patho_afs = np.array(exon_rare_patho_afs)

    pvs1_criteria = np.zeros(len(df), dtype=int)
    # Deal with frameshift/nonsense variants
    non_fs_nmd_variants = (df['Consequence'].str.contains("stop_gained").fillna(False) | \
                           df['Consequence'].str.contains('frameshift').fillna(False)) & \
                           np.logical_not(df['NMD'].fillna(".").str.contains("escaping"))
    logger.info(f"{non_fs_nmd_variants.sum()} variants are frameshift/nonsense that predicted to cause NMD")
    pvs1_criteria[non_fs_nmd_variants & lof_intolerant] = 4

    non_fs_nmd_esp_variants = (df['Consequence'].str.contains("stop_gained").fillna(False) | df['Consequence'].str.contains('frameshift').fillna(False)) & \
                              (df['NMD'].fillna(".").str.contains("escaping") | df['LoF_filter'].fillna(".").str.contains("END_TRUNC"))
    logger.info(f"{non_fs_nmd_esp_variants.sum()} variants are frameshift/nonsense that predicted to escape NMD")
    non_fs_nmd_esp_intol_domain_vars = non_fs_nmd_esp_variants & intolerant_domains
    logger.info(f"{non_fs_nmd_esp_intol_domain_vars.sum()} variants are frameshift/nonsense that predicted to escape NMD and spans functional domains")
    pvs1_criteria[non_fs_nmd_esp_intol_domain_vars & lof_intolerant & (pvs1_criteria < 3)] = 3

    truncate_frac = truncate_fraction(df)
    pvs1_criteria[non_fs_nmd_esp_variants & lof_intolerant & ~intolerant_domains & exon_rare_patho_afs & (truncate_frac >= 0.1) & (pvs1_criteria < 3)] = 3
    pvs1_criteria[non_fs_nmd_esp_variants & lof_intolerant & ~intolerant_domains & exon_rare_patho_afs & (truncate_frac < 0.1) & (pvs1_criteria < 2)] = 2

    # Deal with inframe_deletion variants
    inframe_del_intol_domains = df['Consequence'].str.contains("inframe_deletion") & intolerant_domains
    inframe_ins_intol_domains = df['Consequence'].str.contains("inframe_insertion") & intolerant_domains
    large_inframe_indels = (df["ref"].str.len() > 30) | (df["alt"].str.len() > 30)
    pvs1_criteria[inframe_ins_intol_domains & lof_intolerant & (pvs1_criteria < 3) & ~large_inframe_indels] = 3
    pvs1_criteria[large_inframe_indels & lof_intolerant & (pvs1_criteria < 3) & ~large_inframe_indels] = 3
    pvs1_criteria[df['Consequence'].str.contains("inframe_deletion") & lof_intolerant & ~intolerant_domains & exon_rare_patho_afs & large_inframe_indels & (truncate_frac >= 0.1) & (pvs1_criteria < 3)] = 3
    pvs1_criteria[df['Consequence'].str.contains("inframe_deletion") & lof_intolerant & ~intolerant_domains & exon_rare_patho_afs & large_inframe_indels & (truncate_frac < 0.1) & (pvs1_criteria < 2)] = 2

    # Now we start to deal with splicing variants, we cannot predict NMD escaping given the current annotation, therefore we assume frameshift causing early termination and NMD
    pvs1_criteria[df['splicing_frameshift'] & lof_intolerant & (pvs1_criteria < 4)] = 4
    splicing_inframe_intol_domains = (df['splicing_lof'] | df['splicing_len_changing']) & \
                                      np.logical_not(df['splicing_frameshift'].fillna(False)) & \
                                      df['splicing_span_intol_domain']
    pvs1_criteria[splicing_inframe_intol_domains & lof_intolerant & (pvs1_criteria < 3)] = 3
    pvs1_criteria[(df['splicing_lof'] | df['splicing_len_changing']) & \
                  np.logical_not(df['splicing_frameshift'].fillna(False)) & \
                  np.logical_not(df['splicing_span_intol_domain'].fillna(False)) & \
                  exon_rare_patho_afs & \
                  df['splicing_ten_percent_protein'] & \
                  lof_intolerant & \
                  (pvs1_criteria < 3)] = 3
    pvs1_criteria[(df['splicing_lof'] | df['splicing_len_changing']) & \
                  np.logical_not(df['splicing_frameshift'].fillna(False)) & \
                  np.logical_not(df['splicing_span_intol_domain'].fillna(False)) & \
                  exon_rare_patho_afs & \
                  np.logical_not(df['splicing_ten_percent_protein']) & \
                  clinvar_pathogenic & \
                  (pvs1_criteria < 2)] = 2

    # Thanks to UTRAnnotator, we also implement the PVS1 for frameshift 5UTR variants
    five_utr_frameshift = df['5UTR_frameshift']  #5UTR frameshift bound to NMDs
    five_utr_inframe_intol_domains = df['5UTR_len_changing'] & df['5UTR_span_intol_domain'] & ~five_utr_frameshift
    pvs1_criteria[five_utr_inframe_intol_domains & lof_intolerant & (pvs1_criteria < 3)] = 3
    pvs1_criteria[five_utr_frameshift & lof_intolerant & (pvs1_criteria < 4)] = 4

    alt_start_genes = set(df.groupby(['chrom', 'pos', 'ref', 'alt', 'Gene']).apply(identify_alternative_start_codon_genes).unique().tolist())
    logger.info(f"These are the {len(alt_start_genes)} genes that have functional transcripts using alternative start codons: {alt_start_genes}")
    alt_start_losts = df["Consequence"].str.contains("start_lost") & df["Gene"].isin(alt_start_genes)
    logger.info(f"{alt_start_losts.sum()} variants are having start_lost consequences to transcripts with alternative start codons")
    pvs1_criteria[df['Consequence'].str.contains("start_lost") & ~alt_start_losts & (pvs1_criteria < 2) ] = 2

    return pvs1_criteria, intolerant_domains



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
                        if patho_lp:
                            return False
                        else:
                            return "PS1_Supporting"
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
    elif re.search(r'p\.[A-Z][a-z]+\d+[A-Z][a-z]+$', hgvsp):  # p.Ala123Gly (single amino acid change)
        return "missense"
    else:
        return None


def check_aa_pathogenic(row: dict, 
                        clinvar_tranx_aa_dict: dict, 
                        clinvar_tranx_splice_dict_list: list, 
                        pvs1_strength: int,
                        high_confidence_status = {
                                                    'practice_guideline': 4,                                   # 4 stars
                                                    'reviewed_by_expert_panel': 3,                             # 3 stars
                                                    'criteria_provided,_multiple_submitters,_no_conflicts': 2,  # 2 stars
                                                 }) -> bool:
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
        return check_splice_pathogenic(row, clinvar_tranx_splice_dict_list, pvs1_strength)
        
    # Check if this position has any ClinVar entries
    if raw_protein_pos not in clinvar_tranx_aa_dict:
        logger.debug(f"Protein position {raw_protein_pos} not in ClinVar's VEP annotation records for transcript {transcript}")
        return check_splice_pathogenic(row, clinvar_tranx_splice_dict_list, pvs1_strength)
        
    # Get the clinical significance and review status
    clinvar_entry = clinvar_tranx_aa_dict[raw_protein_pos].get(hgvsp, None)
    
    # Check if any entry is pathogenic with high confidence
    if clinvar_entry:
        logger.debug(f"There is a clinvar entry for {hgvsp} in transcript {transcript} at protein position {raw_protein_pos}, which is {clinvar_entry}")
        for sig, rev_stat in zip(clinvar_entry['CLNSIG'], clinvar_entry['CLNREVSTAT']):
            if (('athogenic' in sig) and (high_confidence_status.get(rev_stat, 0) > 2)) or (("Pathogenic" in sig) and (high_confidence_status.get(rev_stat, 0) >= 2)):
                logger.debug(f"Same_AA_Change: {hgvsp} is pathogenic with high confidence in ClinVar")
                return "Same_AA_Change"
    
    logger.debug(f"No clinvar entry for {hgvsp} in transcript {transcript}. But there are AA changes recorded in the same protein position {raw_protein_pos}")
    splice_pathogenic = check_splice_pathogenic(row, clinvar_tranx_splice_dict_list, pvs1_strength)
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
        
        if query_variant_type != clinvar_variant_type or query_variant_type is None or clinvar_variant_type is None:
            logger.debug(f"Variant types don't match: {hgvsp} ({query_variant_type}) vs {hgvsp_alt} ({clinvar_variant_type})")
            continue
        
        for sig, rev_stat in zip(clinvar_entry['CLNSIG'], clinvar_entry['CLNREVSTAT']):
            if (('athogenic' in sig) and (high_confidence_status.get(rev_stat, 0) > 2)) or (("Pathogenic" in sig) and (high_confidence_status.get(rev_stat, 0) >= 2)):
                logger.debug(f"Same_AA_Residue: {hgvsp} is pathogenic with high confidence in ClinVar and is of the same type as {hgvsp_alt}")
                return "Same_AA_Residue"
            
    return False


def PS1_PM5_criteria(df: pd.DataFrame, 
                     clinvar_aa_dict_pkl: str, 
                     clinvar_splice_dict_pkl: str,
                     ps3_clinvar_patho: np.ndarray,
                     pvs1_criteria: np.ndarray,
                     threads: int = 10, 
                     high_confidence_status = {
                                                'practice_guideline': 4,                                   # 4 stars
                                                'reviewed_by_expert_panel': 3,                             # 3 stars
                                                'criteria_provided,_multiple_submitters,_no_conflicts': 2,  # 2 stars
                                              }) -> Tuple[np.ndarray, np.ndarray]:
    '''
    Identify variants using starmap,
    PS1: Same amino acid change as a previously established pathogenic variant
    PM5: Different amino acid change but same AA residue as a previously established pathogenic variant in a family member
    '''
    logger.info(f"Loading ClinVar AA change dict from {clinvar_aa_dict_pkl}")
    clinvar_aa_dict = pickle.load(gzip.open(clinvar_aa_dict_pkl)) if clinvar_aa_dict_pkl.endswith(".gz") else pickle.load(open(clinvar_aa_dict_pkl, 'rb'))
    logger.info(f"Loading ClinVar splice dict from {clinvar_splice_dict_pkl}")
    clinvar_splice_dict = pickle.load(gzip.open(clinvar_splice_dict_pkl)) if clinvar_splice_dict_pkl.endswith(".gz") else pickle.load(open(clinvar_splice_dict_pkl, 'rb'))
    
    # Convert DataFrame to list of dictionaries
    records = df.to_dict('records')
    
    # Create argument tuples for starmap
    args = [(record, clinvar_aa_dict.get(record['Feature'], {}), clinvar_splice_dict.get(record['Feature'], {}), pvs1_criteria[i], high_confidence_status) for i, record in enumerate(records)]
    
     # Add chunking
    chunk_size = max(len(records) // (threads * 4), 1)
    
    with mp.Pool(threads) as pool:
        logger.info(f"Running check_aa_pathogenic in parallel with {threads} threads on {len(records)} records")
        results = pool.starmap(check_aa_pathogenic, args, chunksize=chunk_size)
    
    results = np.array(results)
    ps1_criteria = (results == "Same_AA_Change") | (results == "PS1")
    ps1_moderate_criteria = (results == "PS1_Moderate")
    ps1_supporting_criteria = (results == "PS1_Supporting")

    pm5_criteria = (results == "Same_AA_Residue") & np.logical_not(df["Consequence"].str.contains("synonymous"))
    pm5_criteria = pm5_criteria & ~ps1_criteria
    ps1_criteria = ps1_criteria & ~ps3_clinvar_patho
    ps1_array = np.zeros(len(df), dtype=int)
    pm5_array = np.zeros(len(df), dtype=int)
    ps1_array[ps1_criteria] = 3
    ps1_array[ps1_moderate_criteria] = 2
    ps1_array[ps1_supporting_criteria] = 1
    pm5_array[pm5_criteria] = 2
    return ps1_array, pm5_array



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
    Currently, we cant implement the variable strength level because the phenotype-gene relationship cant be evaluated by purely bioinformatic means
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
    
    ps2_array = np.zeros(len(df), dtype=int)
    pm6_array = np.zeros(len(df), dtype=int)
    ps2_array[ps2_criteria] = 3
    pm6_array[pm6_criteria] = 2
    
    return ps2_array, pm6_array




def mavedb_interpretation_per_row(row: pd.Series, urn_func_dict=None, mavedb_interpreter=None) -> bool:
    '''
    Interpret the MaveDB scores and return the PS3 and BS3 criteria
    '''
    if pd.isna(row.get('MaveDB_urn', np.nan)):
        row["MaveDB_PS3"] = False
        row["MaveDB_BS3"] = False
        row["MaveDB_score_interpretation"] = np.nan
        return row
    
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


def PS3_BS3_criteria(df: pd.DataFrame, mavedb_metadata_tsv: str = "", high_confidence_status = {
        # Higher confidence (2+ stars)
        'practice_guideline': 4,                                   # 4 stars
        'reviewed_by_expert_panel': 3,                             # 3 stars
        'criteria_provided,_multiple_submitters,_no_conflicts': 2,  # 2 stars
    }) -> pd.DataFrame:
    # Basically rely on ClinVar annotations
    
    result_dict = {"PS3": [], "BS3": [], "clinvar_patho": [], "clinvar_benign": []}

    clinvar_lof = df['CLNSIG'].fillna("").str.contains('Pathogenic') & (df['CLNREVSTAT'].map(high_confidence_status) == 2)
    clinvar_lof = clinvar_lof | (df['CLNSIG'].fillna("").str.contains('athogenic') & (df['CLNREVSTAT'].map(high_confidence_status) >= 3)) # Including Likely_pathogenic
    result_dict["clinvar_patho"] = clinvar_lof

    high_conf_benign = df['CLNSIG'].fillna("").str.contains('Benign') & (df['CLNREVSTAT'].map(high_confidence_status) == 2)
    high_conf_benign = high_conf_benign | (df['CLNSIG'].fillna("").str.contains('enign') & (df['CLNREVSTAT'].map(high_confidence_status) >= 3)) # Including Likely_benign
    result_dict["clinvar_benign"] = high_conf_benign

    ps3_array = np.zeros(len(df), dtype=int)
    bs3_array = np.zeros(len(df), dtype=int)
    
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

    ps3_array[ps3_criteria] = 3
    bs3_array[bs3_criteria] = 3

    result_dict["PS3"] = ps3_array
    result_dict["BS3"] = bs3_array

    return result_dict, df



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




def PP1_criteria(df: pd.DataFrame,
                 recessive: np.ndarray,
                 dominant: np.ndarray,
                 non_monogenic: np.ndarray,
                 non_mendelian: np.ndarray,
                 incomplete_penetrance: np.ndarray,
                 multi_fam_vcf: str = "",
                 multi_fam_ped: str = "",
                 mode: str = "both",) -> np.ndarray:
    '''
    PP1: The variant is cosegregating with a pathogenic variant in one or more families
    '''
    pp1_array = np.zeros(len(df), dtype=int)
    if multi_fam_vcf and multi_fam_ped:
        cosegregating_variants = find_cosegregating_variants(multi_fam_vcf, multi_fam_ped, mode)
        cosegregating_varaints = {mode: {str(t[0]) + ":" + str(t[1]) + ":" + str(t[2]) + "-" + str(t[3]): t[4] for t in variants} for mode,variants in cosegregating_variants.items()}
        recessive_cosegregating = df['variant_id'].map(cosegregating_varaints['recessive'])
        dominant_cosegregating = df['variant_id'].map(cosegregating_varaints['dominant'])
        recessive_ih = np.logical_not(non_monogenic) & np.logical_not(non_mendelian) & np.logical_not(incomplete_penetrance) & recessive
        dominant_ih = np.logical_not(non_monogenic) & np.logical_not(non_mendelian) & np.logical_not(incomplete_penetrance) & dominant & np.logical_not(recessive)

        pp1_array[(recessive_cosegregating <= 2) & recessive_ih] = 1
        pp1_array[(dominant_cosegregating <= 2) & dominant_ih] = 1

        pp1_array[(recessive_cosegregating > 2) & recessive_ih & (recessive_cosegregating <= 4) & (pp1_array < 2)] = 2
        pp1_array[(dominant_cosegregating > 2) & dominant_ih & (dominant_cosegregating <= 4) & (pp1_array < 2)] = 2

        pp1_array[(recessive_cosegregating >= 5) & recessive_ih & (pp1_array < 3)] = 3
        pp1_array[(dominant_cosegregating >= 5) & dominant_ih & (pp1_array < 3)] = 3

        return pp1_array
        
    return pp1_array



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
        


def PM1_criteria(df: pd.DataFrame, 
                 pvs1_criteria: np.ndarray,
                 loc_intol_domain: np.ndarray,
                 intolerant_motifs_pkl: str=None,
                 threads: int = 10) -> np.ndarray:
    # logger.info(f"Loaded the recorded intolerant domains which look alike: {intolerant_domains}")
    row_dicts = df.to_dict('records')
    # args = [(row, intolerant_domains) for row in row_dicts]
    
    # with mp.Pool(threads) as pool:
    #     results = pool.starmap(locate_intolerant_domain, args)

    # loc_intol_domain = np.array(results)
    logger.info(f"There are {np.sum(loc_intol_domain)} variants located in a protein domain that is seemingly intolerant to AA changes according to AM scores")
    
    # We need to check whether the altered splicing event disrupts the intolerant domains downstream of the splicing site, splicing lof already considered affected domains
    splicing_plc_intol_domain = (df['splicing_lof'] | df['splicing_len_changing']) & df['splicing_span_intol_domain']
    utr_plc_intol_domain = (df['5UTR_lof'] | df['5UTR_len_changing']) & df['5UTR_span_intol_domain']
    indel_intol_domain = (df['vep_consq_len_changing'] | df["vep_consq_lof"]) & loc_intol_domain
    nmd_variants = (df["Consequence"].str.contains("stop_gained") | df["Consequence"].str.contains("frameshift")) & np.logical_not(df["NMD"].fillna(".").str.contains("escaping"))
    plc_intol_domain = (splicing_plc_intol_domain | utr_plc_intol_domain | indel_intol_domain) & np.logical_not(nmd_variants)
    indel_no_nmd = (df['vep_consq_len_changing'] | df["vep_consq_lof"]) & np.logical_not(nmd_variants)

    missense = df['Consequence'].str.contains('missense_variant')
    logger.info(f"There are {missense.sum()} missense variants")
    missense_damaging = df["am_class"].fillna("").str.contains('athogenic')
    logger.info(f"There are {missense_damaging.sum()} missense variants that are considered damaging by AlphaMissense")
    missense_damaging = missense_damaging | (df["PrimateAI"] > 0.9)
    logger.info(f"There are {missense_damaging.sum()} missense variants that are considered damaging after considering the PrimateAI score")

    missense_benign = df["am_class"].fillna("").str.contains('benign')
    logger.info(f"There are {missense_benign.sum()} missense variants that are considered benign by AlphaMissense")
    if "MaveDB_PS3" in df.columns:
        # If DMS says the variant is damaging (decreased activity or expression), then consider it as damaging
        missense_damaging = missense_damaging | df["MaveDB_PS3"]
        logger.info(f"There are {missense_damaging.sum()} missense variants that are considered damaging after considering the MaveDB DMS/MPRA data")
    missense_damaging = missense_damaging & missense

    intolerant_motifs = pickle.load(gzip.open(intolerant_motifs_pkl)) if intolerant_motifs_pkl.endswith(".gz") else pickle.load(open(intolerant_motifs_pkl, 'rb'))
    args = [(row, intolerant_motifs.get(row['Feature'], {})) for row in row_dicts]
    logger.info(f"There are {len(args)} variants to be checked for intolerant motifs")
    with mp.Pool(threads) as pool:
        results = pool.starmap(locate_less_char_region, args)
        # Results are a list of tuples, we need to convert them to two independent boolean arrays
        max_am_score_motifs = np.array([r[0] for r in results])
        min_am_score_motifs = np.array([r[1] for r in results])

    logger.info(f"There are {np.sum(min_am_score_motifs)} variants located in a mutational hotspot that is seemingly intolerant to AA changes according to AM scores")
    logger.info(f"There are {np.sum(max_am_score_motifs)} variants located in a mutational hotspot that is seemingly intolerant to AA changes according to most severe AM scores")
    pvs1_double_count = pvs1_criteria >= 3
    pm1_criteria = ( max_am_score_motifs & missense_damaging ) | \
                   ( min_am_score_motifs & missense & np.logical_not(missense_benign) ) | \
                   ( loc_intol_domain & missense & np.logical_not(missense_benign) ) | plc_intol_domain
    pm1_criteria = pm1_criteria & np.logical_not(pvs1_double_count)
    pm1_array = np.zeros(len(df), dtype=int)
    pm1_array[pm1_criteria] = 2
    return pm1_array, loc_intol_domain



def gnomAD_rare_AF(df: pd.DataFrame, cutoff: float) -> np.ndarray:
    if isinstance(cutoff, float):
        return np.where(df['gnomAD_joint_AN_max'].fillna(1000000) >= 1/cutoff, \
                        df['gnomAD_joint_AF_max'].fillna(0) <= cutoff, \
                        df['gnomAD_joint_AF'].fillna(0) <= cutoff)
    elif isinstance(cutoff, str):
        return np.where(df['gnomAD_joint_AN_max'].fillna(1000000) >= 1/(df[cutoff].fillna(10e-6)), \
                        df['gnomAD_joint_AF_max'].fillna(0) <= df[cutoff].fillna(0), \
                        df['gnomAD_joint_AF'].fillna(0) <= df[cutoff].fillna(0))


def gnomAD_rare_nhomalt(df: pd.DataFrame, cutoff: float) -> np.ndarray:
    return np.where(df['gnomAD_joint_AN_max'].fillna(1000000) >= 2/cutoff, \
                    df['gnomAD_nhomalt_max'].fillna(0)/(df['gnomAD_joint_AN_max'].fillna(1000000)/2) <= cutoff, \
                    (df['gnomAD_nhomalt_XX'].fillna(0) + df['gnomAD_nhomalt_XY'].fillna(0))/df['gnomAD_joint_AN'].fillna(1000000) <= cutoff)


def PM2_criteria(df: pd.DataFrame, 
                 clinvar_patho_af_stat: str,
                 clinvar_patho_exon_af_stat: str,
                 gnomAD_extreme_rare_threshold: float = 0.0001) -> np.ndarray:

    logger.info(f"Loading the clinvar pathogenic AF stat from {clinvar_patho_af_stat}")
    clinvar_patho_af_dict = pickle.load(gzip.open(clinvar_patho_af_stat)) if clinvar_patho_af_stat.endswith(".gz") else pickle.load(open(clinvar_patho_af_stat, 'rb'))
    clinvar_patho_af_dict = {k: v for k, v in clinvar_patho_af_dict.items() if v is not None}
    logger.info(f"The clinvar pathogenic AF stat type is {type(clinvar_patho_af_dict)}")
    df.loc[:, "Gene"] = df["Gene"].fillna(np.nan)
    gene_max_patho_af = df['Gene'].map(lambda gene: clinvar_patho_af_dict.get(gene, {"af": 0}).get('af', 0))
    df["clinvar_patho_gene_max_af"] = gene_max_patho_af

    clinvar_patho_exon_af_stat_dict = pickle.load(gzip.open(clinvar_patho_exon_af_stat)) if clinvar_patho_exon_af_stat.endswith(".gz") else pickle.load(open(clinvar_patho_exon_af_stat, 'rb'))
    df["exon_patho_median_af"] = df.apply(lambda row: clinvar_patho_exon_af_stat_dict.get(row['Feature'], {}).get(row["EXON"], (np.nan,))[0], axis=1)
    df["exon_patho_max_af"] = df.apply(lambda row: clinvar_patho_exon_af_stat_dict.get(row['Feature'], {}).get(row["EXON"], (np.nan,np.nan,np.nan))[2], axis=1)

    pm2_supporting_exon = gnomAD_rare_AF(df, "exon_patho_max_af")
    pm2_supporting_gene = gnomAD_rare_AF(df, "clinvar_patho_gene_max_af")

    pm2_supporting = np.where((df["exon_patho_max_af"] == 0) | df["exon_patho_max_af"].isna(), pm2_supporting_gene, pm2_supporting_exon)
    gnomAD_AF_rare = gnomAD_rare_AF(df, gnomAD_extreme_rare_threshold)
    pm2_supporting = np.where((df["clinvar_patho_gene_max_af"] == 0) | df["clinvar_patho_gene_max_af"].isna(), gnomAD_AF_rare, pm2_supporting)

    pm2_moderate = gnomAD_rare_AF(df, "exon_patho_median_af") & gnomAD_AF_rare
    gnomAD_absent = gnomAD_rare_AF(df, 1e-7)
    pm2_moderate = np.where((df["exon_patho_median_af"] == 0) | df["exon_patho_median_af"].isna(), gnomAD_absent, pm2_moderate)

    pm2_array = np.zeros(len(df), dtype=int)
    pm2_array[pm2_supporting] = 1
    pm2_array[pm2_moderate] = 2
    return pm2_array


def PM4_criteria(df: pd.DataFrame, repeat_regions_file: str, loc_intol_domain: np.ndarray) -> np.ndarray:
    # PM4: The variant is causing the protein length change within the Frame
    in_repeat_vars = find_overlaps_bedtools_efficient(df, repeat_regions_file, method="all")
    in_repeat = df['variant_id'].isin(in_repeat_vars)
    frameshift_variants = df["Consequence"].str.contains("frameshift")
    nmd_variants = df["Consequence"].str.contains("stop_gained") & np.logical_not(df["NMD"].fillna(".").str.contains("escaping"))
    inframe_indels = (df['vep_consq_len_changing'] & \
                     (np.logical_not(in_repeat) | loc_intol_domain) & \
                     np.logical_not(frameshift_variants) & \
                     np.logical_not(nmd_variants))
    splicing_inframe_indels = df["splicing_len_changing"] & \
                              np.logical_not(df["splicing_frameshift"]) & \
                              (np.logical_not(in_repeat) | df["splicing_span_intol_domain"])
    utr_inframe_indels = df["5UTR_len_changing"] & \
                         (np.logical_not(in_repeat) | df["5UTR_span_intol_domain"])
    pm4_criteria = inframe_indels | splicing_inframe_indels | utr_inframe_indels
    pm4_array = np.zeros(len(df), dtype=int)
    pm4_array[pm4_criteria] = 2
    return pm4_array, in_repeat



def analyze_bp1_pp2(gene_stat_dict: Dict):
    """
    Analyze BP1 and PP2 criteria for a gene based on variant statistics.
    
    Args:
        gene_stat_dict: Dictionary containing gene-level variant statistics
        
    Returns:
        Dictionary with BP1 and PP2 analysis results
    """
    gene = gene_stat_dict['ensembl_id']
    if len(gene_stat_dict) <= 1:
        return (gene, (False, False))

    # Get pathogenic variants
    patho_variants = gene_stat_dict['clinvar_pathogenic']
    total_patho = len(patho_variants)
    
    if total_patho == 0:
        return (gene, (False, False))
    
    # BP1 Analysis
    # Combine large AA change and NMD variants
    large_aa_nmd_variants = (gene_stat_dict['large_aachange'] | 
                            gene_stat_dict['putative_nmd_variants'])
    patho_large_aa_nmd = len(patho_variants & large_aa_nmd_variants)
    bp1_fraction = patho_large_aa_nmd / total_patho
    bp1_granted = bp1_fraction >= 0.7
    
    # PP2 Analysis
    # Create contingency table for small vs large AA changes
    small_aa_variants = gene_stat_dict['small_aachange']
    large_aa_variants = gene_stat_dict['large_aachange']
    
    # Count variants in each category
    patho_small_aa = len(patho_variants & small_aa_variants)
    patho_large_aa = len(patho_variants & large_aa_variants)
    non_patho_small_aa = len(small_aa_variants - patho_variants)
    non_patho_large_aa = len(large_aa_variants - patho_variants)
    
    # Create contingency table
    table = [[patho_small_aa, patho_large_aa],
             [non_patho_small_aa, non_patho_large_aa]]
    
    # Initialize PP2 results
    pp2_granted = False
    pp2_test_used = 'none'
    pp2_pvalue = np.nan
    pp2_odds_ratio = np.nan
    
    # Check if we can perform statistical tests
    total_count = sum(sum(row) for row in table)
    min_cell_count = min(min(row) for row in table)
    
    if (total_count >= 5 and 
        all(sum(row) > 0 for row in table) and 
        all(sum(col) > 0 for col in zip(*table))):
        
        if total_count >= 20 and min_cell_count >= 5:
            # Use Chi-square test
            chi2, pp2_pvalue = stats.chi2_contingency(table)[0:2]
            pp2_odds_ratio = ((table[0][1] * table[1][0]) / 
                            (table[0][0] * table[1][1]))  # Note: reversed for large AA changes
            pp2_test_used = 'chi_square'
        else:
            # Use Fisher's exact test
            pp2_odds_ratio, pp2_pvalue = stats.fisher_exact(table, alternative='greater')
            pp2_test_used = 'fisher'
        
        # Adjust p-value of 1 to 0.99999 to avoid log10(1) = 0
        pp2_pvalue = min(pp2_pvalue, 0.99999)
        
        # Grant PP2 if there is NO significant enrichment of large AA changes
        # (p-value > 0.05) and missense variants are ≥50% of pathogenic variants
        pp2_fraction = patho_small_aa / total_patho
        pp2_granted = (pp2_pvalue > 0.05) and (pp2_fraction >= 0.4)
    else:
        # Fall back to fraction-based approach
        pp2_fraction = patho_small_aa / total_patho
        pp2_granted = pp2_fraction >= 0.5
        pp2_test_used = 'fraction'
    
    result_dict = {
        'bp1_granted': bp1_granted,
        'bp1_fraction': bp1_fraction,
        'pp2_granted': pp2_granted,
        'pp2_fraction': patho_small_aa / total_patho if total_patho > 0 else 0.0,
        'pp2_test_used': pp2_test_used,
        'pp2_pvalue': pp2_pvalue,
        'pp2_odds_ratio': pp2_odds_ratio
    }

    logger.debug(f"The result_dict looks like {result_dict}")

    return (gene, (result_dict['pp2_granted'], result_dict['bp1_granted']))



def PP2_BP1_criteria(df: pd.DataFrame, 
                     clinvar_stats_dict: dict,
                     am_intol_domains_tsv: str,
                     threads: int = 10) -> Tuple[pd.Series, pd.Series]:
    """
    Efficient implementation of PP2/BP1 criteria evaluation using 
    ClinVar statistics dictionary.
    
    Args:
        df: Variant annotation DataFrame
        clinvar_stats_dict: Dictionary containing ClinVar statistics per gene
        threads: Number of threads for multiprocessing
        
    Returns:
        Tuple[pd.Series, pd.Series]: PP2 and BP1 criteria boolean series
    """
    
    # Get unique genes
    unique_genes = df['Gene'].dropna().unique()
    logger.info(f"Processing BP1/PP2 criteria for {len(unique_genes)} unique genes")
    
    # Process genes in parallel
    with mp.Pool(threads) as pool:
        gene_results = pool.imap_unordered(analyze_bp1_pp2, (clinvar_stats_dict.get(gene, {"ensembl_id": gene}) for gene in unique_genes))
        # Convert to dictionary for faster lookups
        gene_criteria_dict = dict(gene_results)
    
    # Create result Series using map operation (vectorized)
    gene_pp2 = df['Gene'].map(lambda g: gene_criteria_dict.get(g, (False, False))[0] if pd.notna(g) else False)
    gene_bp1 = df['Gene'].map(lambda g: gene_criteria_dict.get(g, (False, False))[1] if pd.notna(g) else False)
    
    # Apply PP2 only to missense variants
    is_missense = df['Consequence'].str.contains('missense_variant', na=False)
    pp2_criteria = gene_pp2 & is_missense
    
    # BP1 applies to all variants
    bp1_criteria = gene_bp1 & is_missense
    
    # Log stats
    logger.info(f"Granted PP2 to {pp2_criteria.sum()} variants based on ClinVar stats")
    logger.info(f"Granted BP1 to {bp1_criteria.sum()} variants based on ClinVar stats")

    # Load the AlphaMissense intolerant domains
    am_intol_domains_df = pd.read_table(am_intol_domains_tsv, low_memory=False)
    am_intol_domains_df["Ensembl_Gene_ID"] = am_intol_domains_df["domain"].str.split(":").str[0]
    by_gene = am_intol_domains_df.groupby("Ensembl_Gene_ID", as_index=False)
    gene_all_intol_domains = {}
    for gene, gene_df in by_gene:
        if gene_df["is_more_tolerant"].any():
            gene_all_intol_domains[gene] = False
        else:
            gene_all_intol_domains[gene] = True

    # Apply PP2 to all variants
    pp2_am = df['Gene'].map(lambda g: gene_all_intol_domains.get(g, False))
    pp2_gene_am = df["Gene_avg_AM_score"].fillna(0) > 0.564
    pp2_criteria = pp2_criteria | pp2_am | pp2_gene_am
    pp2_criteria = pp2_criteria & is_missense
    logger.info(f"Granted PP2 to {pp2_criteria.sum()} variants based on updates from AlphaMissense intolerant domains")

    pp2_array = np.zeros(len(df), dtype=int)
    pp2_array[pp2_criteria] = 1
    bp1_array = np.zeros(len(df), dtype=int)
    bp1_array[bp1_criteria] = 1
    return pp2_array, bp1_array



def PP3_BP4_criteria(df: pd.DataFrame, high_confidence_status = {
        # Higher confidence (2+ stars)
        'practice_guideline': 4,                                   # 4 stars
        'reviewed_by_expert_panel': 3,                             # 3 stars
        'criteria_provided,_multiple_submitters,_no_conflicts': 2,  # 2 stars
    }) -> pd.Series:
    # PP3: predicted to be deleterious by in-silico tools
    # Including PrimateAI (Missense), CADD, AlphaMissense (Missense), VEP, SpliceAI, SpliceVault, UTRAnnotator
    
    missense_variant = df['Consequence'].str.contains('missense_variant') & np.logical_not(df['Consequence'].str.contains('stop')) & np.logical_not(df['Consequence'].str.contains('frameshift'))
    splice_variant = df['Consequence'].str.contains('splice') & np.logical_not(df['Consequence'].str.contains('stop')) & np.logical_not(df['Consequence'].str.contains('frameshift'))
    five_utr_variant = df['Consequence'].str.contains('5_prime_UTR_variant').fillna(False)
    missense_variant = missense_variant.fillna(False)
    splice_variant = splice_variant.fillna(False)
    # BP4: variant is reported benign
    pp3_criteria = ((df['PrimateAI'].fillna(0).astype(float) > 0.9) & missense_variant) | \
                    ((df['CADD_phred'].fillna(10).astype(float) >= 20) & np.logical_not(splice_variant) & np.logical_not(missense_variant) & np.logical_not(five_utr_variant)) | \
                    (df['am_class'].fillna("").str.contains('pathogenic') & missense_variant) | \
                    (df['vep_consq_lof'] & np.logical_not(splice_variant) & np.logical_not(missense_variant) & np.logical_not(five_utr_variant)) | \
                    (((df['splicing_lof'] | (df['CADD_phred'].fillna(10).astype(float) >= 20)) & splice_variant) | (df['5UTR_lof'] & five_utr_variant))
    clinvar_benign = df['CLNSIG'].fillna("").str.contains('enign') & (df['CLNREVSTAT'].map(high_confidence_status, na_action="ignore") >= 2)
    pp3_criteria = pp3_criteria & ~clinvar_benign
    
    missense_benign = (df['PrimateAI'].fillna(0).astype(float) < 0.8).fillna(True) & (df['am_pathogenicity'].fillna(0).astype(float) < 0.564).fillna(True) & missense_variant
    splice_benign = np.logical_not(df['splicing_lof'].fillna(False)) & splice_variant & (df['CADD_phred'].fillna(10).astype(float) < 20)
    utr_benign = np.logical_not(df['5UTR_lof'].fillna(False)) & five_utr_variant
    other_benign = np.logical_not(df['vep_consq_lof'].fillna(False)) & (df['CADD_phred'].fillna(10).astype(float) < 20).fillna(True) & (df["CADD_reg_phred"].fillna(10).astype(float) < 20).fillna(True) & np.logical_not(splice_variant) & np.logical_not(missense_variant) & np.logical_not(five_utr_variant)
    bp4_criteria = missense_benign | splice_benign | utr_benign | other_benign
    clinvar_patho = df['CLNSIG'].fillna("").str.contains('athogenic') & (df['CLNREVSTAT'].map(high_confidence_status, na_action="ignore") >= 2)
    bp4_criteria = bp4_criteria & ~clinvar_patho
    
    pp3_array = np.zeros(len(df), dtype=int)
    bp4_array = np.zeros(len(df), dtype=int)
    pp3_array[pp3_criteria] = 1
    bp4_array[bp4_criteria] = 1
    return pp3_array, bp4_array
           
   

def PP5_BP6_criteria(df: pd.DataFrame, clinvar_patho, clinvar_benign) -> pd.Series:
    # PP5: The variant is reported as pathogenic by a reputable source but without to many supporting evidences
    pp5_criteria = df['CLNSIG'].fillna("").str.contains('athogenic') & ~clinvar_patho
    bp6_criteria = df['CLNSIG'].fillna("").str.contains('enign') & ~clinvar_benign

    pp5_array = np.zeros(len(df), dtype=int)
    bp6_array = np.zeros(len(df), dtype=int)
    pp5_array[pp5_criteria] = 1
    bp6_array[bp6_criteria] = 1
    return pp5_array, bp6_array



def BS1_criteria(df: pd.DataFrame, 
                 gene_to_am_score_map: dict,
                 expected_incidence: float = 0.001,
                 gene_dosage_sensitivity: str = "",
                 pm2_criteria: np.ndarray = None,
                 threads: int = 10):
    # BS1: PAF of variant is greater than expected incidence of the disease
    autosomal = (df['chrom'] != "chrX") & (df['chrom'] != "chrY")
    x_linked = df['chrom'] == "chrX"
    y_linked = df['chrom'] == "chrY"

    recessive, dominant, non_monogenic, non_mendelian, haplo_insufficient, incomplete_penetrance = identify_inheritance_mode(df, gene_to_am_score_map, gene_dosage_sensitivity, threads)

    false_neg_rate, common_vars = control_false_neg_rate(df['gnomAD_joint_AF_max'], df['gnomAD_joint_AN_max'], af_threshold=expected_incidence, alpha=0.01)

    max_ind_incidence = np.where(df['gnomAD_joint_AN_max']/2 > 10/expected_incidence, df['gnomAD_nhomalt_max']/(df['gnomAD_joint_AN_max']/2), (df['gnomAD_nhomalt_XX'] + df['gnomAD_nhomalt_XY'])/(df['gnomAD_joint_AN']/2))
    max_af_larger_incidence = np.where(common_vars.isna(), df['gnomAD_joint_AF'] > expected_incidence, common_vars)
    logger.info(f"There are {max_af_larger_incidence.sum()} variants having their PAF greater than the expected incidence of the disease")
    # For autosomal dominant disease, we can assign BS1 if the variant is observed the frequency of the variant is greater than the expected incidence of the disease
    autosomal_dominant = autosomal & dominant & max_af_larger_incidence & np.logical_not(recessive) & np.logical_not(non_monogenic) & np.logical_not(non_mendelian) & np.logical_not(incomplete_penetrance)
    autosomal_recessive = autosomal & recessive & (max_ind_incidence > expected_incidence) & np.logical_not(non_monogenic) & np.logical_not(non_mendelian) & np.logical_not(incomplete_penetrance)

    # For X-linked disease, we can assign BS1 if the variant is observed the frequency of the variant is greater than the expected incidence of the disease
    x_linked_recessive = x_linked & recessive & ((df['gnomAD_nhomalt_XX'] + df['gnomAD_nhomalt_XY'])/(df['gnomAD_joint_AN']/2) > expected_incidence) & np.logical_not(non_monogenic) & np.logical_not(non_mendelian) & np.logical_not(incomplete_penetrance)
    x_linked_dominant = x_linked & dominant & max_af_larger_incidence & np.logical_not(non_monogenic) & np.logical_not(non_mendelian) & np.logical_not(incomplete_penetrance)
    # For Y-linked disease, we can assign BS1 if the variant is observed the frequency of the variant is greater than the expected incidence of the disease
    y_linked = y_linked & max_af_larger_incidence & np.logical_not(non_monogenic) & np.logical_not(non_mendelian) & np.logical_not(incomplete_penetrance)
    greater_than_disease_incidence = autosomal_dominant | autosomal_recessive | x_linked_recessive | x_linked_dominant | y_linked
    logger.info(f"There are {greater_than_disease_incidence.sum()} variants having their PAF greater than the expected incidence of the disease")
    gene_max_patho_af = df["clinvar_patho_gene_max_af"].fillna(0)

    _, greater_than_clinvar_patho_af = control_false_neg_rate(df['gnomAD_joint_AF_max'], df['gnomAD_joint_AN_max'], af_threshold=gene_max_patho_af, alpha=0.01)
    greater_than_clinvar_patho_af = np.where(greater_than_clinvar_patho_af.isna(), df['gnomAD_joint_AF'] > gene_max_patho_af, greater_than_clinvar_patho_af)
    _, greater_than_basic_af = control_false_neg_rate(df['gnomAD_joint_AF_max'], df['gnomAD_joint_AN_max'], af_threshold=0.0001, alpha=0.01)
    greater_than_basic_af = np.where(greater_than_basic_af.isna(), df['gnomAD_joint_AF'] > 0.0001, greater_than_basic_af)
    bs1_criteria = (greater_than_disease_incidence | greater_than_clinvar_patho_af) & greater_than_basic_af
    bs1_criteria = bs1_criteria & (pm2_criteria == 0)
    bs1_array = np.zeros(len(df), dtype=int)
    bs1_array[bs1_criteria] = 3
    return bs1_array, recessive, dominant, non_monogenic, non_mendelian, haplo_insufficient, incomplete_penetrance



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
    recessive_set = {"Autosomal recessive inheritance", "X-linked recessive inheritance"}

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


def identify_inheritance_mode_per_row(row_dict: dict, gene_mean_am_score: float, clingen_curate_score: int = None) -> Tuple[bool, bool]:
    # We need to use three fields of the table to determine the inheritance mode:
    # 1. Gene
    # 2. LOEUF
    # 3. HPO_IDs
    # 4. HPO_gene_inheritance (overrides the above two fields), HPO observed dominant inheritance can derive from GOF variants
    # 5. ClinGen curated dosage sensitivity, 3 means haploinsufficient, 30 or 40 means haplosufficient

    loeuf_score = float(row_dict.get('LOEUF', 0.3))
    loeuf_score = 1.0 if pd.isna(loeuf_score) else loeuf_score  # If LOEUF is NaN, we leave the decision to gene avg AM score
    haplo_insufficient = (loeuf_score <= 0.36) or (gene_mean_am_score >= 0.564)
    haplo_insufficient = haplo_insufficient and (loeuf_score <= 0.8) and (gene_mean_am_score >= 0.4)
    haplo_sufficient = not haplo_insufficient

    hpo_inheritance = parse_hpo_inheritance(row_dict)
    if isinstance(hpo_inheritance, bool):
        logger.debug(f"No HPO inheritance information for {row_dict['Gene']}, using LOEUF: {row_dict['LOEUF']} and AM score: {gene_mean_am_score} to determine inheritance mode. The haploinsufficiency is {haplo_insufficient} and haplosufficiency is {haplo_sufficient}")
        return haplo_sufficient, haplo_insufficient, False, False, haplo_insufficient, hpo_inheritance

    if clingen_curate_score:
        logger.debug(f"Using ClinGen curated dosage sensitivity to determine inheritance mode for {row_dict['Gene']}, the ClinGen record looks like this: \n{clingen_curate_score}\n")
        if clingen_curate_score == 3:
            hpo_inheritance['hpo_recessive'] = False
            hpo_inheritance['hpo_dominant'] = True
        elif clingen_curate_score == 30 or clingen_curate_score == 40:
            hpo_inheritance['hpo_recessive'] = True
            # Cannot modify hpo_dominant because it might relates to GOF variants 

    if hpo_inheritance['hpo_recessive']:
        haplo_insufficient = False

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



def BS2_criteria(df: pd.DataFrame, 
                 recessive: np.ndarray,
                 dominant: np.ndarray,
                 non_monogenic: np.ndarray,
                 non_mendelian: np.ndarray,
                 incomplete_penetrance: np.ndarray,
                 pm2_criteria: np.ndarray) -> pd.Series:
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
    autosomal = (df['chrom'] != "chrX") & (df['chrom'] != "chrY") & (df['chrom'] != "chrM")
    x_linked = df['chrom'] == "chrX"
    y_linked = df['chrom'] == "chrY"

    # HP:0003829: Incomplete penetrance
    # HP:4000159: Moderate penetrance
    # HP:4000160: Low penetrance
    not_late_onsets = df['HPO_IDs'].map(hpo_onset_modes)
    complete_penetrance = df['HPO_IDs'].fillna("").str.contains("HP:0034950") # HP:0034950: Complete penetrance
    moderate_to_low_penetrance = df['HPO_IDs'].fillna("").str.contains("HP:4000159") | df['HPO_IDs'].fillna("").str.contains("HP:4000160")
    
    # For autosomal dominant disease, we can assign BS2 if the variant is observed in a healthy adult (either homozygous or heterozygous)
    autosomal_dominant = autosomal & dominant & ((df['gnomAD_joint_AF'].fillna(0) * df['gnomAD_joint_AN'].fillna(0)) > 5) & np.logical_not(recessive) & np.logical_not(non_monogenic) & np.logical_not(non_mendelian) & np.logical_not(moderate_to_low_penetrance)
    autosomal_recessive = autosomal & recessive & (((df['gnomAD_nhomalt_XX'].fillna(0) + df['gnomAD_nhomalt_XY'].fillna(0)) > 5)) & np.logical_not(non_monogenic) & np.logical_not(non_mendelian) & np.logical_not(moderate_to_low_penetrance)

    # For X-linked disease, we can assign BS2 if the variant is observed in a healthy adult male (hemizygous) or a healthy adult female (homozygous)
    x_linked_recessive = x_linked & recessive & (((df['gnomAD_nhomalt_XX'].fillna(0) + df['gnomAD_nhomalt_XY'].fillna(0)) > 5)) & np.logical_not(non_monogenic) & np.logical_not(non_mendelian) & np.logical_not(moderate_to_low_penetrance)
    x_linked_dominant = x_linked & dominant & np.logical_not(recessive) & (((df['gnomAD_nhomalt_XX'].fillna(0) + df['gnomAD_nhomalt_XY'].fillna(0)) > 5) | ((df['gnomAD_joint_AF'].fillna(0) * df['gnomAD_joint_AN'].fillna(0)) > 5)) & np.logical_not(non_monogenic) & np.logical_not(non_mendelian) & np.logical_not(moderate_to_low_penetrance)

    y_linked = y_linked & (df['gnomAD_joint_AF_XY'].fillna(0) > 0) & np.logical_not(non_monogenic) & np.logical_not(non_mendelian) & np.logical_not(moderate_to_low_penetrance)
    
    if 'control_AC' in df.columns:
        logger.warning("Seems user has provided a control cohort VCF, using control allele counts to assign BS2 criteria")
        autosomal_dominant = (autosomal_dominant | (df['control_AC'].fillna(0) > 0)) & np.logical_not(moderate_to_low_penetrance)
        x_linked_dominant = (x_linked_dominant | (df['control_AC'].fillna(0) > 0)) & np.logical_not(moderate_to_low_penetrance)
        y_linked = (y_linked | (df['control_AC'].fillna(0) > 0)) & np.logical_not(moderate_to_low_penetrance)
    if 'control_nhomalt' in df.columns:
        logger.warning("Seems user has provided a control cohort VCF, using control homozygous counts to assign BS2 criteria")
        autosomal_recessive = (autosomal_recessive | (df['control_nhomalt'].fillna(0) > 0)) & np.logical_not(moderate_to_low_penetrance)
        x_linked_recessive = (x_linked_recessive | (df['control_nhomalt'].fillna(0) > 0)) & np.logical_not(moderate_to_low_penetrance)

    gnomad_bs2_criteria = autosomal_dominant | autosomal_recessive | x_linked_recessive | x_linked_dominant | y_linked
    bs2_criteria = gnomad_bs2_criteria & np.logical_not(incomplete_penetrance) & not_late_onsets

    # For combination of complete penetrance and not late onset, any AC or nhomalt in gnomAD is considered BS2
    autosomal_dominant = autosomal & dominant & ((df['gnomAD_joint_AF'].fillna(0) * df['gnomAD_joint_AN'].fillna(0)) > 0) & np.logical_not(recessive) & np.logical_not(non_monogenic) & np.logical_not(non_mendelian)
    autosomal_recessive = autosomal & recessive & (((df['gnomAD_nhomalt_XX'].fillna(0) + df['gnomAD_nhomalt_XY'].fillna(0)) > 0)) & np.logical_not(non_monogenic) & np.logical_not(non_mendelian)
    x_linked_recessive = x_linked & recessive & (((df['gnomAD_nhomalt_XX'].fillna(0) + df['gnomAD_nhomalt_XY'].fillna(0)) > 0)) & np.logical_not(non_monogenic) & np.logical_not(non_mendelian)
    x_linked_dominant = x_linked & dominant & np.logical_not(recessive) & (((df['gnomAD_nhomalt_XX'].fillna(0) + df['gnomAD_nhomalt_XY'].fillna(0)) > 0) | ((df['gnomAD_joint_AF'].fillna(0) * df['gnomAD_joint_AN'].fillna(0)) > 0)) & np.logical_not(non_monogenic) & np.logical_not(non_mendelian)

    bs2_complete_penetrance = (autosomal_dominant | autosomal_recessive | x_linked_recessive | x_linked_dominant | y_linked) & complete_penetrance & not_late_onsets
    bs2_criteria = bs2_complete_penetrance | bs2_criteria
    bs2_criteria = bs2_criteria & (pm2_criteria == 0)
    bs2_array = np.zeros(len(df), dtype=int)
    bs2_array[bs2_criteria] = 3
    return bs2_array




def BS4_criteria(df: pd.DataFrame, ped_df: pd.DataFrame, fam_name: str, 
                 recessive: np.ndarray, 
                 dominant: np.ndarray, 
                 non_monogenic: np.ndarray, 
                 non_mendelian: np.ndarray, 
                 incomplete_penetrance: np.ndarray) -> pd.Series:
    # BS4: lack of family segregation
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
        hmem_criteria = ((df[healthy_mem].str.count("1") == 2) & recessive) | \
                        ((df[healthy_mem].str.count("1") >= 1) & dominant & np.logical_not(recessive))
        hmem_criteria = hmem_criteria & np.logical_not(non_monogenic) & np.logical_not(non_mendelian) & np.logical_not(incomplete_penetrance)
        final_criteria = final_criteria | hmem_criteria

    bs4_array = np.zeros(len(df), dtype=int)
    bs4_array[final_criteria] = 3
    return bs4_array


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
                     pm2_criteria: np.ndarray,
                     pvs1_criteria: np.ndarray,
                     ps1_criteria: np.ndarray,
                     ps3_criteria: np.ndarray,
                     is_recessive: np.ndarray,
                     is_dominant: np.ndarray,
                     incomplete_penetrance: np.ndarray) -> Tuple[pd.Series, pd.Series]:
    # BP2: observed in trans with a pathogenic variant in dominant disease, Or in-cis with a pathogenic variant with any inheritance mode
    # PM3: observed in trans with a pathogenic variant in recessive disease.
    pathogenic = df["vep_consq_lof"] | df["splicing_lof"] | df["5UTR_lof"] | (ps1_criteria & ps3_criteria) | pvs1_criteria

    in_cis_pathogenic, in_trans_pathogenic = batch_annotate_cis_trans_from_table(df, 
                                                                                 pathogenic, 
                                                                                 ped_df, 
                                                                                 chrom_col="chrom",
                                                                                 pos_col="pos",
                                                                                 ref_col="ref",
                                                                                 alt_col="alt",
                                                                                 gene_col="Gene",
                                                                                 ped_sample_id_col="IndividualID", 
                                                                                 ped_paternal_id_col="PaternalID", 
                                                                                 ped_maternal_id_col="MaternalID", 
                                                                                 ped_sex_col="Sex", 
                                                                                 ped_phenotype_col="Phenotype", 
                                                                                 ped_patient_value=2)

    in_cis_pathogenic = in_cis_pathogenic > 0
    in_trans_pathogenic = in_trans_pathogenic > 0
    
    bp2_criteria = (in_trans_pathogenic & is_dominant & np.logical_not(is_recessive)) | (in_cis_pathogenic & is_recessive)
    pm3_criteria = in_trans_pathogenic & (is_recessive) & pm2_criteria

    bp2_criteria = bp2_criteria & ~incomplete_penetrance
    pm3_criteria = pm3_criteria & ~incomplete_penetrance

    bp2_array = np.zeros(len(df), dtype=int)
    pm3_array = np.zeros(len(df), dtype=int)
    bp2_array[bp2_criteria] = 1
    pm3_array[pm3_criteria] = 2

    return bp2_array, pm3_array


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


def BP3_criteria(df: pd.DataFrame, repeat_region_file: str, interpro_entry_map_pkl: str, pm1_criteria: np.ndarray, pm4_criteria: np.ndarray) -> pd.Series:
    # BP3: in-frame deletion in a repetitive region without a known function
    inframe_del = (df['Consequence'].fillna("").str.contains('inframe_deletion')) | \
                  (df['Consequence'].fillna("").str.contains('inframe_insertion')) | \
                  (df['vep_consq_len_changing'] & ~df["Consequence"].str.contains("frameshift")) | \
                  (df['splicing_len_changing'] & ~df['splicing_frameshift']) | \
                  (df['5UTR_len_changing'] & ~df['5UTR_frameshift']) | \
                  (df["Consequence"].fillna("").str.contains("stop_gained") & df["NMD"].fillna("").str.contains("escaping"))

    # Repeat region file is a gzipped bed file, we can read it with pandas
    in_repeat_regions = find_overlaps_bedtools_efficient(df, repeat_region_file)
    dm_instance = DomainNormalizer()
    interpro_entry_map_dict = pickle.load(gzip.open(interpro_entry_map_pkl)) if interpro_entry_map_pkl.endswith(".gz") else pickle.load(open(interpro_entry_map_pkl, "rb"))
    repetitive_region = ( df['DOMAINS'].fillna(".").str.contains('Low_complexity') | df['variant_id'].isin(in_repeat_regions) )
    repetitive_region = repetitive_region & np.logical_not(df["5UTR_span_intol_domain"]) & np.logical_not(df["splicing_span_intol_domain"])
    bp3_criteria = inframe_del & repetitive_region & (pm1_criteria == 0) & (pm4_criteria == 0)
    bp3_array = np.zeros(len(df), dtype=int)
    bp3_array[bp3_criteria] = 1
    return bp3_array


def BP5_criteria(df: pd.DataFrame, 
                 alt_disease_vcf: str, 
                 gene_to_am_score_map: dict, 
                 clingen_dosage_sensitivity: str,
                 threads: int = 10) -> np.ndarray:
    '''
    BP5: variant found in a sample with known alternative molecular basis for disease
    Considers inheritance mode when determining presence in alternative disease samples
    '''
    # Convert DataFrame to dictionaries for parallel processing
    row_dicts = df.to_dict('records')

    clingen_dosage_df = pd.read_table(clingen_dosage_sensitivity, low_memory=False).dropna(subset=["#Gene Symbol", "Haploinsufficiency Score"])
    clingen_dosage_map = dict(zip(clingen_dosage_df['#Gene Symbol'], clingen_dosage_df['Haploinsufficiency Score'].astype(int)))
    
    # Prepare arguments for starmap
    args = [(row_dict, alt_disease_vcf, gene_to_am_score_map.get(row_dict['Gene'], np.nan), clingen_dosage_map.get(row_dict['SYMBOL'], None)) for row_dict in row_dicts]
    
    # Process in parallel
    with mp.Pool(threads) as pool:
        results = pool.starmap(identify_presence_in_alt_disease_vcf, args)
    
    bp5_criteria = np.array(results)
    bp5_array = np.zeros(len(df), dtype=int)
    bp5_array[bp5_criteria] = 1
    return bp5_array
    

def identify_presence_in_alt_disease_vcf(row: dict, alt_disease_vcf: str, gene_mean_am_score: float, clingen_curate_score: int) -> bool:
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
        is_recessive, is_dominant, is_non_monogenic, is_non_mendelian, haplo_insufficient, incomplete_penetrance = identify_inheritance_mode_per_row(row, gene_mean_am_score, clingen_curate_score)
        
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
    bp7_array = np.zeros(len(df), dtype=int)
    bp7_array[synonymous & no_splicing_altering & not_conserved] = 1
    return bp7_array


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


def create_criteria_summary(row, criteria_order=None, clingen_evidence=None):
    """
    Create a summary string of active criteria for a given row.
    
    Args:
        row: A row from the criteria matrix
        criteria_order: List of criteria in order
        
    Returns:
        A semicolon-separated string of active criteria
    """

    strength_level_suffix = {1: "Supporting", 2: "Moderate", 3: "Strong", 4: "Very Strong", 5: "Standalone"}
    label_strength_level = {"P": "Supporting", "M": "Moderate", "S": "Strong", "V": "Very Strong", "A": "Standalone"}

    active_criteria = [ col if strength_level_suffix.get(int(row[col]), None) == label_strength_level.get(col[1], None) else col + "_" + strength_level_suffix.get(int(row[col]), None) for col in criteria_order if int(row[col]) > 0 ]
    clingen_criteria = clingen_evidence.get(row["variant_id"], None)
    if clingen_criteria:
        logger.info(f"ClinGen evidence for variant {row['variant_id']}: {clingen_criteria}")
        active_criteria = [ fix_clingen_term(c) for c in clingen_criteria.split(",") ] if isinstance(clingen_criteria, str) else []

    return ";".join(active_criteria) if active_criteria else np.nan


def fix_clingen_term(clingen_str: str) -> str:
    if not isinstance(clingen_str, str):
        raise ValueError(f"ClinGen evidence input {clingen_str} is not a string")
    
    if "_" not in clingen_str:
        return clingen_str
    else:
        prefix = clingen_str.split("_")[0]
        suffix = clingen_str.split("_")[1]
        if suffix.startswith("Very"):
            suffix = "Very Strong"
        elif suffix.startswith("Stand"):
            suffix = "Stand Alone"
        return f"{prefix}_{suffix}"


def translate_strength_level(criteria_name: str, criteria_list_str: str) -> int:
    strength_map = { "V": 4, "S": 3, "M": 2, "P": 1, "A": 5}
    variable_strength_level = { "Supporting": 1, "Moderate": 2, "Strong": 3, "Very Strong": 4, "Stand Alone": 5, "Very": 4, "Stand": 5, "Standalone": 5 }
    if not isinstance(criteria_list_str, str):
        return 0
    elif criteria_name not in criteria_list_str:
        return 0
    else:
        criteria_list = criteria_list_str.split(";")
        if criteria_name in criteria_list:
            return strength_map[criteria_name[1]]
        else:
            criteria_match = [c for c in criteria_list if c.startswith(criteria_name)][0]
            return variable_strength_level[criteria_match.split("_")[1]]


def summarize_acmg_criteria(df: pd.DataFrame, criteria_dict: Dict[str, np.ndarray], clingen_map_pkl: str) -> Tuple[pd.DataFrame, pd.DataFrame]:
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
        {name: criteria_dict.get(name, np.zeros(len(df), dtype=int)) 
         for name in criteria_order},
        index=df.index
    )

    criteria_matrix["variant_id"] = df["variant_id"]
    logger.info(f"Criteria matrix: \n{criteria_matrix[:20].to_string(index=False)}")

    # Check which column has no integer values and return error
    for col in criteria_order:
        if not criteria_matrix[col].isin([0, 1, 2, 3, 4, 5]).all():
            raise ValueError(f"Column {col} has non-integer values")

    # Load ClinGen evidence map
    clingen_map = pickle.load(gzip.open(clingen_map_pkl, "rb")) if clingen_map_pkl.endswith(".gz") else pickle.load(open(clingen_map_pkl, "rb"))
    
    # Add summary column to original DataFrame
    df['ACMG_criteria'] = criteria_matrix.apply(create_criteria_summary, axis=1, criteria_order=criteria_order, clingen_evidence=clingen_map)

    # We need to use the ACMG_criteria value to translate to the criteria_matrix
    clingen_match = criteria_matrix["variant_id"].isin(clingen_map)
    for col in criteria_order:
        criteria_matrix.loc[clingen_match, col] = df.loc[clingen_match, "ACMG_criteria"].map(lambda x: translate_strength_level(col, x))
    
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

    # evidence_weights = {"PVS1": 1,
    #                     "PS1": 1/exp_base, "PS2": 1/exp_base, "PS3": 1/exp_base, "PS4": 1/exp_base,
    #                     "PM1": 1/exp_base**2, "PM2": 1/exp_base**3, "PM3": 1/exp_base**2, "PM4": 1/exp_base**2, "PM5": 1/exp_base**2, "PM6": 1/exp_base**2,
    #                     "PP1": 1/exp_base**3, "PP2": 1/exp_base**3, "PP3": 1/exp_base**3, "PP4": 1/exp_base**3, "PP5": 1/exp_base**3,
    #                     "BA1": -5,
    #                     "BS1": -1/exp_base, "BS2": -1/exp_base, "BS3": -1/exp_base, "BS4": -1/exp_base,
    #                     "BP1": -1/exp_base**3, "BP2": -1/exp_base**3, "BP3": -1/exp_base**3, "BP4": -1/exp_base**3, "BP5": -1/exp_base**3, "BP6": -1/exp_base**3, "BP7": -1/exp_base**3}

    evidence_weights = {4: 1.0, 3: 1/exp_base, 2: 1/exp_base**2, 1: 1/exp_base**3, 5: 2.0}

    # Default odds function (Illustrative - replace with a calibrated function based on SVI guidelines)
    odds_function = lambda total_weight: np.power(odds_pvst, total_weight)

    # There are 2 internal inconsistencies in the SVI guidelines.
    # Pathogenic(ii) and Likely Pathogenic(i) do not generate a posterior probability of the corresponding class.
    # We need to handle this by setting the odds of pathogenic to 0 if the total weight is negative. If we run into 2 PS evidences, we adjust the total_weight_pathogenic to 1514
    # If there is a PVS + 1 PM, we need to adjust their sum to 350

    patho_matches = [n for c,n in criteria_assignment.items() if c.startswith("P") and int(n) > 0]
    ps_matches = [n for n in patho_matches if n == 3]
    pm_matches = [n for n in patho_matches if n == 2]
    pp_matches = [n for n in patho_matches if n == 1]
    pvs_matches = [n for n in patho_matches if n == 4]
    benign_matches = [n for c,n in criteria_assignment.items() if c.startswith("B") and int(n) > 0]
    bs_matches = [n for n in benign_matches if n == 3]
    bp_matches = [n for n in benign_matches if n == 1]
    ba_matches = [n for n in benign_matches if n == 5]
    bm_matches = [n for n in benign_matches if n == 2]

    odds_benign = odds_function(-float(sum([evidence_weights[n] for n in bs_matches + bp_matches + ba_matches + bm_matches])))
    if len(ps_matches) >= 2:
        ps_odds = 1514 * odds_function(float(sum([evidence_weights[int(n)] for n in ps_matches[:-2]])))
    else:
        ps_odds = odds_function(float(sum([evidence_weights[int(n)] for n in ps_matches])))
    
    if len(pvs_matches) > 0 and len(pm_matches) > 0 and len(pp_matches) == 0:
        pv_pm_odds = odds_pvst * odds_function(float(sum([evidence_weights[int(n)] for n in pm_matches[:-1]])))
    else:
        pv_pm_odds = odds_function(float(sum([evidence_weights[int(n)] for n in pvs_matches + pm_matches])))

    pp_odds = odds_function(float(sum([evidence_weights[int(n)] for n in pp_matches])))

    odds_path = ps_odds * pp_odds * pv_pm_odds

    # Calculate odds of pathogenicity
    odds_path = odds_path * odds_benign

    # Calculate posterior probability
    posterior_probability = (odds_path * prior_probability) / ((odds_path - 1) * prior_probability + 1)
    # Round posterior probability to 3 decimal places
    posterior_probability = round(posterior_probability, 4)

    if len(ba_matches) > 0:
        posterior_probability = 0

    # Classify the variant based on the posterior probability
    if posterior_probability > 0.99:
        acmg_class = "Pathogenic"
    elif posterior_probability > 0.899:
        acmg_class = "Likely Pathogenic"
    elif posterior_probability < 0.0005:
        acmg_class = "Benign"
    elif posterior_probability < 0.05:
        acmg_class = "Likely Benign"
    else:
        acmg_class = "Uncertain Significance"

    return posterior_probability, acmg_class




def sort_and_rank_variants(df: pd.DataFrame, ped_df: pd.DataFrame, fam_name: str, pvs1_criteria: np.ndarray, gene_to_am_score_map: dict, recessive: np.ndarray, dominant: np.ndarray) -> pd.DataFrame:
    """
    Sort variants by their maximum ACMG quantitative score and add ranking.
    
    Args:
        df: DataFrame with ACMG_quant_score column and variant information
        
    Returns:
        DataFrame sorted by variant's max ACMG score with added variant rank column
    """
    if ped_df is not None and fam_name is not None:
        proband = ped_df.loc[(ped_df['#FamilyID'] == fam_name) & (ped_df['Phenotype'].isin(["2", 2])), 'IndividualID'].values[0]
        proband_het = (df.loc[:, proband].str.count("1") == 1)
    else:
        proband_het = np.array([False] * len(df))
        
    df["sort_index"] = df.loc[:, "ACMG_quant_score"]

    # normalized_loeuf = df.loc[:, "LOEUF"].fillna(1)
    # normalized_mean_am = (1-df.loc[:, "Gene"].map(gene_to_am_score_map)) * 2
    # normalized_mean_am = normalized_mean_am.fillna(1)
    haplo_sufficient = (df.loc[:, "LOEUF"] > 0.8) | (df.loc[:, "Gene"].map(gene_to_am_score_map) < 0.4)
    df["haplo_insuf_index"] = 1
    df.loc[haplo_sufficient, "haplo_insuf_index"] = 0.9 # This penalty coefficient will lead to a near-1 posterior prob downgrade to below 0.9 which is below the cutoff of likely pathogenic.
    
    only_recessive = (recessive & np.logical_not(dominant))
    lof = pvs1_criteria | df["vep_consq_lof"] | df["splicing_lof"] | df["5UTR_lof"]
    df.loc[lof & proband_het & only_recessive, "sort_index"] = df.loc[lof & proband_het & only_recessive, "sort_index"] * df.loc[lof & proband_het & only_recessive, "haplo_insuf_index"]
    # Group by variant coordinates to get max score per variant
    variant_groups = df.groupby(['chrom', 'pos', 'ref', 'alt'])
    max_scores = variant_groups['sort_index'].transform('max')
    
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

    ba1_array = np.zeros(len(anno_df), dtype=int)
    ba1_array[ba1_criteria] = 5
    return ba1_array



def ACMG_criteria_assign(anno_table: str, 
                         am_score_table: str, 
                         clinvar_patho_af_stat: str,
                         clinvar_patho_exon_af_stat: str,
                         clinvar_aa_dict_pkl: str,
                         clinvar_splice_dict_pkl: str,
                         interpro_entry_map_pkl: str,
                         intolerant_domains_pkl: str,
                         am_intol_domains_tsv: str,
                         intolerant_motifs_pkl: str,
                         clinvar_gene_stat_pkl: str,
                         tranx_exon_domain_map_pkl: str,
                         repeat_region_file: str,
                         gene_dosage_sensitivity: str,
                         pp1_vcf: str = "",
                         pp1_ped: str = "",
                         mavedb_metadata_tsv: str = "",
                         fam_name: str = "",
                         ped_table: str = "",
                         alt_disease_vcf: str = "",
                         clingen_map_pkl: str = "",
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
    2. PP4, Patient's phenotype or family history is highly specific for a disease with a single genetic etiology. (This cannot be applied simutaneously with PP1, because if a disease is tightly linked with only one gene, the segregation is doomed. Therefore further segregation observation does not add any more confidence because the confidence from this perspective is already reaching a ceiling.)

    ===========Regarding the criteria, Can be applied===========
    1. PVS1, LOEUF <= 0.35 should be considered as intolerant to LoF.
    2. PM2, is reduced from Moderate to Supporting.
    3. PP5/BP6, reputable source reported as benign or pathogenic. Suggested to be abandoned or at least not assigned along with PS3/BS3.
    3. PM3, can be only applied if PM2 is True (sufficiently rare in gnomAD).
    """
    anno_df = pd.read_table(anno_table, low_memory=False).drop_duplicates()
    logger.info(f"Got {threads} threads to process the input table {anno_table}, now the table looks like: \n{anno_df[:5].to_string(index=False)}")
    anno_df = vep_consq_interpret(anno_df, threads)

    am_score_df = pd.read_table(am_score_table, low_memory=False)
    ped_df = None
    fam_df = None
    proband = None
    if ped_table:
        ped_df = pd.read_table(ped_table, low_memory=False)
        fam_df = ped_df[ped_df['#FamilyID'] == fam_name]
        proband = fam_df.loc[fam_df['Phenotype'].isin(["2", 2]), 'IndividualID'].values[0]
        

    # Convert the am_score_df to a dictionary:
    # 1. Ensembl transcript ID (column 'transcript') to mean AM score (column 'mean_am_pathogenicity')
    am_score_dict = dict(zip(am_score_df['transcript'], am_score_df['mean_am_pathogenicity']))
    # Use anno_df to create a dict map from Ensembl transcript ID to gene ID
    transcript_to_gene_map = dict(zip(anno_df['Feature'], anno_df['Gene']))
    # Use the two dict above to create dict that maps gene ID to mean AM score
    gene_to_am_score_map = {g: am_score_dict[t] for t, g in transcript_to_gene_map.items() if t in am_score_dict}
    clinvar_aa_dict = pickle.load(gzip.open(clinvar_aa_dict_pkl)) if clinvar_aa_dict_pkl.endswith(".gz") else pickle.load(open(clinvar_aa_dict_pkl, "rb"))
    clinvar_aa_gene_map = {g: clinvar_aa_dict[t] for t, g in transcript_to_gene_map.items() if t in clinvar_aa_dict}
    logger.info(f"gene_to_am_score_map created, {len(gene_to_am_score_map)} genes are having the AM score")

    # Establish the variant ID column
    anno_df["variant_id"] = anno_df["chrom"] + ":" + anno_df["pos"].astype(str) + ":" + anno_df["ref"] + "-" + anno_df["alt"]

    # Load the intolerant domains
    if intolerant_domains_pkl.endswith(".gz"):
        with gzip.open(intolerant_domains_pkl, 'rb') as f:
            mm = mmap.mmap(f.fileno(), 0, access=mmap.ACCESS_READ)
            intolerant_domains = pickle.load(mm)
    else:
        with open(intolerant_domains_pkl, 'rb') as f:
            mm = mmap.mmap(f.fileno(), 0, access=mmap.ACCESS_READ)
            intolerant_domains = pickle.load(mm)
    logger.info(f"Loading the recorded intolerant domains which look alike: {intolerant_domains}")

    # Apply the PVS1 criteria, LoF on a gene known to to be pathogenic due to LoF
    pvs1_criteria, locate_intol_domains = PVS1_criteria(anno_df, 
                                                        clinvar_aa_gene_map, 
                                                        clinvar_patho_exon_af_stat,
                                                        interpro_entry_map_pkl,
                                                        gene_to_am_score_map,
                                                        intolerant_domains=intolerant_domains,
                                                        tranx_exon_domain_map_pkl=tranx_exon_domain_map_pkl,
                                                        proband_gt_col=proband ) # When test on ClinVar variants, fam_name is set to None because no genotype information are provided
    anno_df["span_functional_domains"] = locate_intol_domains
    logger.info(f"PVS1 criteria applied, {(pvs1_criteria > 0).sum()} variants are having the PVS1 criteria")
    gc.collect()

    # Apply the PS3 and BS3 criteria
    ps3bs3_results, anno_df = PS3_BS3_criteria(anno_df, mavedb_metadata_tsv)
    ps3_criteria = ps3bs3_results['PS3']
    bs3_criteria = ps3bs3_results['BS3']
    logger.info(f"PS3 criteria applied, {(ps3_criteria > 0).sum()} variants are having the PS3 criteria")
    logger.info(f"BS3 criteria applied, {(bs3_criteria > 0).sum()} variants are having the BS3 criteria")
    gc.collect()

    # Apply the PS1 and PM5 criteria
    ps1_criteria, pm5_criteria = PS1_PM5_criteria(anno_df, clinvar_aa_dict_pkl, clinvar_splice_dict_pkl, ps3bs3_results['clinvar_patho'], pvs1_criteria, threads)
    logger.info(f"PS1 criteria applied, {(ps1_criteria > 0).sum()} variants are having the PS1 criteria") 
    logger.info(f"PM5 criteria applied, {(pm5_criteria > 0).sum()} variants are having the PM5 criteria")
    gc.collect()

    # Apply the PS2 and PM6 criteria
    if not ped_df is None and not fam_name is None:
        ps2_criteria, pm6_criteria = PS2_PM6_criteria(anno_df, ped_df, fam_name, threads=threads)
    else:
        logger.warning(f"No ped_table provided, skip the PS2 and PM6 criteria")
        ps2_criteria, pm6_criteria = np.array([0] * len(anno_df)), np.array([0] * len(anno_df))
    logger.info(f"PS2 criteria applied, {(ps2_criteria > 0).sum()} variants are having the PS2 criteria")
    logger.info(f"PM6 criteria applied, {(pm6_criteria > 0).sum()} variants are having the PM6 criteria")
    gc.collect()
   
    '''
    PS4 cannot be applied because usually we dont have enough cases to determine the frequency of the variant
    '''
    # Apply PM1 criteria, mutational hotspot or well-established functional protein domain
    pm1_criteria, loc_intol_domain = PM1_criteria(anno_df, pvs1_criteria, locate_intol_domains, intolerant_motifs_pkl, threads)
    # pm1_criteria = pm1_criteria & ~ps1_criteria  # PS1 is already a strength to indicate the intolerance to the AA changes incurred by the variant
    logger.info(f"PM1 criteria applied, {(pm1_criteria > 0).sum()} variants are having the PM1 criteria")
    gc.collect()
    # Apply PM4 criteria, causing the protein length change
    pm4_criteria, in_repeat_vars = PM4_criteria(anno_df, repeat_region_file, loc_intol_domain)
    anno_df["variant_in_repeat_region"] = in_repeat_vars
    gc.collect()
    logger.info(f"PM4 criteria applied, {(pm4_criteria > 0).sum()} variants are having the PM4 criteria")

    # Apply PP2 criteria, missense variant in a gene/domain that not only intolerant to truncating variants but also intolerant to missense variants
    clinvar_gene_stat = pickle.load(gzip.open(clinvar_gene_stat_pkl, "rb")) if clinvar_gene_stat_pkl.endswith(".gz") else pickle.load(open(clinvar_gene_stat_pkl, "rb"))
    pp2_criteria, bp1_criteria = PP2_BP1_criteria(anno_df, clinvar_gene_stat, am_intol_domains_tsv)
    logger.info(f"PP2 criteria applied, {(pp2_criteria > 0).sum()} variants are having the PP2 criteria")
    logger.info(f"BP1 criteria applied, {(bp1_criteria > 0).sum()} variants are having the BP1 criteria")
    gc.collect()
    # Apply PP3 criteria, predicted to be deleterious by in-silico tools
    pp3_criteria, bp4_criteria = PP3_BP4_criteria(anno_df)
    # bp4_criteria = bp4_criteria & ~bs3_criteria
    logger.info(f"BP4 criteria applied, {(bp4_criteria > 0).sum()} variants are having the BP4 criteria")
    logger.info(f"PP3 criteria applied, {(pp3_criteria > 0).sum()} variants are having the PP3 criteria")
    gc.collect()
    '''
    PP4 cannot be applied, Patient's phenotype or family history is highly specific for a disease with a single genetic etiology
    '''
    # Apply PP5 criteria, reported as pathogenic by a reputable source but without to many supporting evidences
    # Apply BP6 criteria, reported as benign by a reputable source but without to many supporting evidences
    pp5_criteria, bp6_criteria = PP5_BP6_criteria(anno_df, ps3bs3_results['clinvar_patho'], ps3bs3_results['clinvar_benign'])
    logger.info(f"PP5 criteria applied, {(pp5_criteria > 0).sum()} variants are having the PP5 criteria")
    logger.info(f"BP6 criteria applied, {(bp6_criteria > 0).sum()} variants are having the BP6 criteria")
    gc.collect()
    
    ba1_criteria = BA1_criteria(anno_df)

    # Apply PM2 criteria, absent from gnomAD or extremely rare in gnomAD
    pm2_criteria = PM2_criteria(anno_df, 
                                clinvar_patho_af_stat,
                                clinvar_patho_exon_af_stat,
                                gnomAD_extreme_rare_threshold)
    logger.info(f"PM2 criteria applied, {(pm2_criteria > 0).sum()} variants are having the PM2 criteria")
    gc.collect()

    # Apply BS1, PAF of variant is greater than expected incidence of the disease
    bs1_criteria, recessive, dominant, non_monogenic, non_mendelian, haplo_insufficient, incomplete_penetrance = BS1_criteria(anno_df, 
                                                                                                                              gene_to_am_score_map, 
                                                                                                                              threads = threads, 
                                                                                                                              expected_incidence = expected_incidence, 
                                                                                                                              gene_dosage_sensitivity = gene_dosage_sensitivity,
                                                                                                                              pm2_criteria = pm2_criteria)
    logger.info(f"BS1 criteria applied, {(bs1_criteria > 0).sum()} variants are having the BS1 criteria")
    gc.collect()

    # Summarize the inheritance mode, first prepare a df composed of recessive, dominant, non_monogenic, non_mendelian, haplo_insufficient, incomplete_penetrance
    inheritance_df = pd.DataFrame({
        'recessive': recessive,
        'dominant': dominant,
        'non_monogenic': non_monogenic,
        'non_mendelian': non_mendelian,
        'haplo_insufficient': haplo_insufficient,
        'incomplete_penetrance': incomplete_penetrance
    })

    inheritance_df.loc[recessive, "recessive"] = "recessive"
    inheritance_df.loc[~recessive, "recessive"] = ""
    inheritance_df.loc[dominant, "dominant"] = "dominant"
    inheritance_df.loc[~dominant, "dominant"] = ""
    inheritance_df.loc[non_monogenic, "non_monogenic"] = "non_monogenic"
    inheritance_df.loc[~non_monogenic, "non_monogenic"] = ""
    inheritance_df.loc[non_mendelian, "non_mendelian"] = "non_mendelian"
    inheritance_df.loc[~non_mendelian, "non_mendelian"] = ""
    inheritance_df.loc[incomplete_penetrance, "incomplete_penetrance"] = "incomplete_penetrance"
    inheritance_df.loc[~incomplete_penetrance, "incomplete_penetrance"] = ""
    inheritance_df.loc[haplo_insufficient, "haplo_insufficient"] = "haplo_insufficient"
    inheritance_df.loc[~haplo_insufficient, "haplo_insufficient"] = ""
    anno_df["inheritance_mode"] = inheritance_df.apply(lambda x: ",".join([i for i in x if i != ""]), axis=1)

    # Apply BS2, variant observed in a healthy adult
    bs2_criteria = BS2_criteria(anno_df, recessive, dominant, non_monogenic, non_mendelian, incomplete_penetrance, pm2_criteria)
    logger.info(f"BS2 criteria applied, {(bs2_criteria > 0).sum()} variants are having the BS2 criteria")
    gc.collect()

    # Apply PP1 criteria, variant is cosegregating with a pathogenic variant in one or more families
    pp1_criteria = PP1_criteria(anno_df, recessive, dominant, non_monogenic, non_mendelian, incomplete_penetrance, pp1_vcf, pp1_ped)
    logger.info(f"PP1 criteria applied, {(pp1_criteria > 0).sum()} variants are having the PP1 criteria")
    gc.collect()

    # Apply BS4, lack of family segregation
    if not ped_df is None and not fam_name is None:
        bs4_criteria = BS4_criteria(anno_df, ped_df, fam_name, recessive, dominant, non_monogenic, non_mendelian, incomplete_penetrance)
    else:
        logger.warning(f"No ped_table provided, skip the BS4 criteria")
        bs4_criteria = np.array([0] * len(anno_df))
    logger.info(f"BS4 criteria applied, {(bs4_criteria > 0).sum()} variants are having the BS4 criteria")
    # Apply BP3, in-frame indels in a repetitive region without a known function
    bp3_criteria = BP3_criteria(anno_df, repeat_region_file, interpro_entry_map_pkl, pm1_criteria, pm4_criteria)
    logger.info(f"BP3 criteria applied, {(bp3_criteria > 0).sum()} variants are having the BP3 criteria")
    gc.collect()
    # Apply BP5, variant found in a sample with known alternative molecular basis for disease
    if alt_disease_vcf:
        bp5_criteria = BP5_criteria(anno_df, alt_disease_vcf, gene_to_am_score_map, gene_dosage_sensitivity, threads)
    else:
        bp5_criteria = np.array([0] * len(anno_df))
    logger.info(f"BP5 criteria applied, {(bp5_criteria > 0).sum()} variants are having the BP5 criteria")
    gc.collect()
    # Apply BP7, synonymous variant, no splicing-altering consequence, not conserved. 
    bp7_criteria = BP7_criteria(anno_df)
    logger.info(f"BP7 criteria applied, {(bp7_criteria > 0).sum()} variants are having the BP7 criteria")
    gc.collect()

    # Apply BP2, observed in trans with a pathogenic variant in dominant disease, Or in-cis with a pathogenic variant with any inheritance mode
    # Apply PM3, observed in trans with a pathogenic variant in recessive disease.
    if not ped_df is None and not fam_name is None:
        bp2_criteria, pm3_criteria = BP2_PM3_criteria(anno_df, 
                                                      ped_df, 
                                                      pm2_criteria,
                                                      pvs1_criteria,
                                                      ps1_criteria,
                                                      ps3_criteria,
                                                      recessive,
                                                      dominant,
                                                      incomplete_penetrance)
    else:
        logger.warning(f"No ped_table provided, skip the BP2 and PM3 criteria")
        bp2_criteria, pm3_criteria = np.array([0] * len(anno_df)), np.array([0] * len(anno_df))
    logger.info(f"BP2 criteria applied, {(bp2_criteria > 0).sum()} variants are having the BP2 criteria")
    logger.info(f"PM3 criteria applied, {(pm3_criteria > 0).sum()} variants are having the PM3 criteria")
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
        'PP1': pp1_criteria,
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
    anno_df, criteria_matrix = summarize_acmg_criteria(anno_df, criteria_dict, clingen_map_pkl)
    # Save the criteria matrix to a file which the path is based on the input anno_table
    output_matrix = ".".join(anno_table.split(".")[:-1]) + ".acmg.tsv"
    criteria_matrix.to_csv(output_matrix, sep="\t", index=False)
    
    # Apply quantification using ACMG_criteria column and use that to sort the variants
    posterior_probability, acmg_class = zip(*criteria_matrix.apply(calculate_posterior_probability, axis=1))
    anno_df["ACMG_quant_score"] = posterior_probability
    anno_df["ACMG_class"] = acmg_class

    # Sort and rank variants
    anno_df = sort_and_rank_variants(anno_df, ped_df, fam_name, pvs1_criteria, gene_to_am_score_map, recessive, dominant)
    # Save the annotated table to replace the input anno_table
    anno_df.to_csv(anno_table, sep="\t", index=False)
    
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
    parser.add_argument("--gene_dosage_sensitivity", type=str, required=True)
    parser.add_argument("--am_intol_domains_tsv", type=str, required=True)
    parser.add_argument("--intolerant_motifs_pkl", type=str, required=True)
    parser.add_argument("--clinvar_gene_stat_pkl", type=str, required=True)
    parser.add_argument("--tranx_exon_domain_map_pkl", type=str, required=True)
    parser.add_argument("--am_score_vcf", type=str, required=False, default=None)
    parser.add_argument("--alt_disease_vcf", type=str, required=False, default=None)
    parser.add_argument("--clingen_map_pkl", type=str, required=False, default=None)
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
                                                    args.am_intol_domains_tsv,
                                                    args.intolerant_motifs_pkl,
                                                    args.clinvar_gene_stat_pkl,
                                                    args.tranx_exon_domain_map_pkl,
                                                    args.repeat_region_file,
                                                    args.gene_dosage_sensitivity,
                                                    mavedb_metadata_tsv=args.mavedb_metadata_tsv,
                                                    fam_name=args.fam_name,
                                                    ped_table=args.ped_table,
                                                    alt_disease_vcf=args.alt_disease_vcf,
                                                    clingen_map_pkl=args.clingen_map_pkl,
                                                    gnomAD_extreme_rare_threshold=args.gnomAD_extreme_rare_threshold,
                                                    expected_incidence=args.expected_incidence,
                                                    threads=args.threads)




