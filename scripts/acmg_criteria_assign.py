#!/usr/bin/env python

import pysam
import pickle
import logging
import pandas as pd
import numpy as np
import argparse as ap
from collections import namedtuple
from typing import Tuple, Dict
import multiprocessing as mp
from multiprocessing import Manager

from stat_protein_domain_amscores import nested_defaultdict


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
    
    is_lof = any(c in consq for c in lof_criteria)
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



def PVS1_criteria(df: pd.DataFrame, am_score_dict: dict) -> pd.DataFrame:
    # LoF is high confidence if it is a LoF variant and the consequence is known
    # VEP LoF, splicing LoF, UTRAnnotator LoF, ClinVar Pathogenic
    high_confidence_status = {
        # Higher confidence (2+ stars)
        'practice_guideline',                                    # 4 stars
        'reviewed_by_expert_panel',                             # 3 stars
        'criteria_provided,_multiple_submitters,_no_conflicts',  # 2 stars
    }

    clinvar_lof = (df['CLNSIG'] == 'Pathogenic') & df['CLNREVSTAT'].isin(high_confidence_status)
    lof_criteria = df['vep_consq_lof'] | df['splicing_lof'] | df['5UTR_lof'] | clinvar_lof

    # Determine whether the gene is known to be pathogenic due to LoF
    # Evaluated by mean AM scores and LOEUF
    am_intolerant_tranx = df['Feature'].map(am_score_dict).fillna(0) >= 0.58  # The bin 0 threshold according to AM publication
    loeuf_intolerant_tranx = df['LOEUF'] < 0.6 # The most two intolerant bins according to LOEUF publication
    intolerant_lof = am_intolerant_tranx | loeuf_intolerant_tranx
    pvs1_criteria = lof_criteria & intolerant_lof

    return pvs1_criteria



def check_aa_pathogenic(row: dict, clinvar_tranx_aa_dict: dict, high_confidence_status: set) -> bool:
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
    protein_pos = str(raw_protein_pos).split("/")[0] if raw_protein_pos and not raw_protein_pos in [np.nan, np.inf] else '' # e.g "117"
    hgvsp = row.get('HGVSp', '') # e.g ENSP00000349098.5:p.E117K
    
    logger.debug(f"The current row records a variant overlapping with transcript {transcript} at protein position {protein_pos} with HGVSp {hgvsp}")
    if hgvsp in [np.nan, np.inf, 'nan', 'inf', '']:
        return False
        
    # Check if this transcript has any ClinVar entries
    if not clinvar_tranx_aa_dict:
        logger.debug(f"Transcript {transcript} not in ClinVar's VEP annotation records")
        return False
        
    # Check if this position has any ClinVar entries
    if protein_pos not in clinvar_tranx_aa_dict:
        logger.debug(f"Protein position {protein_pos} not in ClinVar's VEP annotation records for transcript {transcript}")
        return False
        
    # Get the clinical significance and review status
    clinvar_entry = clinvar_tranx_aa_dict[protein_pos].get(hgvsp, None)
    
    # Check if any entry is pathogenic with high confidence
    if clinvar_entry:
        logger.info(f"There is a clinvar entry for {hgvsp} in transcript {transcript} at protein position {protein_pos}")
        for sig, rev_stat in zip(clinvar_entry['CLNSIG'], clinvar_entry['CLNREVSTAT']):
            if ('athogenic' in sig and  # Captures both 'Pathogenic' and 'Likely_pathogenic'
                any(status in rev_stat for status in high_confidence_status)):
                logger.info(f"Same_AA_Change: {hgvsp} is pathogenic with high confidence in ClinVar")
                return "Same_AA_Change"
    
    logger.debug(f"No clinvar entry for {hgvsp} in transcript {transcript}. But there are AA changes recorded in the same protein position {protein_pos}")
    for hgvsp_alt, clinvar_entry in clinvar_tranx_aa_dict[protein_pos].items():
        for sig, rev_stat in zip(clinvar_entry['CLNSIG'], clinvar_entry['CLNREVSTAT']):
            if ('athogenic' in sig and  # Captures both 'Pathogenic' and 'Likely_pathogenic'
                any(status in rev_stat for status in high_confidence_status)):
                logger.info(f"Same_AA_Residue: {hgvsp} is pathogenic with high confidence in ClinVar")
                return "Same_AA_Residue"
            
    return False


def PS1_PM5_criteria(df: pd.DataFrame, 
                     clinvar_aa_dict_pkl: str, 
                     threads: int = 10) -> Tuple[np.ndarray, np.ndarray]:
    '''
    Identify variants using starmap,
    PS1: Same amino acid change as a previously established pathogenic variant
    PM5: Different amino acid change but same AA residue as a previously established pathogenic variant in a family member
    '''
    high_confidence_status = {
        'practice_guideline',
        'reviewed_by_expert_panel',
        'criteria_provided,_multiple_submitters,_no_conflicts'
    }
    clinvar_aa_dict = pickle.load(open(clinvar_aa_dict_pkl, 'rb'))
    
    # Convert DataFrame to list of dictionaries
    records = df.to_dict('records')
    
    # Create argument tuples for starmap
    args = [(record, clinvar_aa_dict.get(record['Feature'], {}), high_confidence_status) for record in records]
    
     # Add chunking
    chunk_size = max(len(records) // (threads * 4), 1)
    
    with mp.Pool(threads) as pool:
        logger.info(f"Running check_aa_pathogenic in parallel with {threads} threads on {len(records)} records")
        results = pool.starmap(check_aa_pathogenic, args, chunksize=chunk_size)
    
    results = np.array(results)
    ps1_criteria = results == "Same_AA_Change"
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
            elif row.get('gnomAD_joint_af') in [0, np.nan]:
                return "PM6"
            else:
                return False
    else:
        if row['chrom'] == "chrX":
            if "1" in row.get(father, '') or "1" in row.get(mother, ''):
                return False
            elif father and mother:
                return "PS2"
            elif row.get('gnomAD_joint_af') in [0, np.nan]:
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
            elif row.get('gnomAD_joint_af') in [0, np.nan]:
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
    proband_info, father_info, mother_info, sib_info = identify_fam_members(ped_df, fam_name)
    
    # Convert DataFrame to list of dictionaries for multiprocessing
    records = df.to_dict('records')
    
    # Prepare arguments for starmap
    args = [(record, (proband_info, father_info, mother_info, sib_info), ped_df.loc[ped_df['#FamilyID'] == fam_name, :]) for record in records]
    
    # Use pool.starmap with the arguments
    with mp.Pool(threads) as pool:
        results = pool.starmap(determine_denovo_per_row, args)
    
    results = np.array(results)
    ps2_criteria = results == "PS2"
    pm6_criteria = results == "PM6"
    
    return ps2_criteria, pm6_criteria




def PS3_BS3_criteria(df: pd.DataFrame) -> pd.DataFrame:
    # Basically rely on ClinVar annotations
    high_confidence_status = {
        # Higher confidence (2+ stars)
        'practice_guideline',                                    # 4 stars
        'reviewed_by_expert_panel',                             # 3 stars
        'criteria_provided,_multiple_submitters,_no_conflicts',  # 2 stars
    }

    clinvar_lof = df['CLNSIG'].str.contains('athogenic') & df['CLNREVSTAT'].isin(high_confidence_status)
    high_conf_benign = df['CLNSIG'].str.contains('enign') & df['CLNREVSTAT'].isin(high_confidence_status)

    return {
            'PS3': clinvar_lof,
            'BS3': high_conf_benign
            }


def locate_intolerant_domain(row: dict, intolerant_domains: set) -> bool:
    '''
    Check if the variant is located in a mutational hotspot or a well-established functional protein domain
    The intolerant domains are calculated by scripts: 
    1. scripts/am_pick_intolerant_domains.py
    2. scripts/clinvar_pick_intolerant_domains.py
    '''
    gene = row.get('Gene', '') if isinstance(row.get('Gene', None), str) else ''
    if isinstance(row.get('DOMAINS', None), str):
        domains = [gene + "_" + d.replace(":", "_") for d in row.get('DOMAINS', None).split('&')]
    else:
        return False
    # The variant itself need to be protein altering.
    missense = row['Consequence'] == 'missense_variant'
    length_changing = row['vep_consq_length_changing'] | row['splicing_len_changing'] | row['5UTR_len_changing']
    protein_altering = missense | length_changing | (row['Consequence'] == 'protein_altering_variant')

    return any(domain in intolerant_domains for domain in domains) & protein_altering
    


def PM1_criteria(df: pd.DataFrame, 
                 intolerant_domains_pkl: str, 
                 threads: int = 10) -> np.ndarray:
    # PM1: The variant is located in a mutational hotspot or a well-established functional protein domain
    # Load intolerant domains into shared memory
    row_dicts = df.to_dict('records')
    with Manager() as manager:
        shared_domains = manager.dict()
        shared_domains['intolerant_domains'] = pickle.load(open(intolerant_domains_pkl, 'rb'))
        
        args = [(row, shared_domains['intolerant_domains']) for row in row_dicts]
        
        with mp.Pool(threads) as pool:
            results = pool.starmap(locate_intolerant_domain, args)
    
    return np.array(results)


def PM2_criteria(df: pd.DataFrame, gnomAD_extreme_rare_threshold: float = 0.0001) -> np.ndarray:
    # PM2: The variant is absent from gnomAD or the variant is extremely rare in gnomAD
    gnomAD_absent = (df['gnomAD_joint_af'] == 0) | (df['gnomAD_joint_af'].isna())
    gnomAD_rare = df['gnomAD_joint_af'] < gnomAD_extreme_rare_threshold
    return gnomAD_absent | gnomAD_rare


def PM4_criteria(df: pd.DataFrame) -> np.ndarray:
    # PM4: The variant is causing the protein length change
    return df['vep_consq_length_changing'] | df['splicing_len_changing'] | df['5UTR_len_changing']



def PP2_BP1_criteria(df: pd.DataFrame, domain_mechanism_tsv: str) -> Tuple[pd.Series, pd.Series]:
    '''
    Evaluate PP2 criteria: missense variant in a gene that has a low rate of benign missense variants
    and where missense variants are a common mechanism of disease.

    Simutaenously evaluating BP1 criteria: missense variant located in a gene/domain where primarily intolerant to truncating variants 
    
    Args:
        df: Variant annotation DataFrame
        domain_mechanism_tsv: Path to domain mechanism analysis results
        
    Returns:
        Boolean Series indicating PP2 criteria met
    '''
    # Read domain mechanism data
    domain_df = pd.read_table(domain_mechanism_tsv, low_memory=False)
    
    # Filter for genes and domains that are both intolerant to truncating and missense variants
    tolerant_genes = set(domain_df.loc[domain_df['gene_pvalue'] > 0.05,'gene_id'].unique())
    tolerant_domains = set(domain_df.loc[domain_df['pvalue'] > 0.05, 'domain'].unique())
    intolerant_genes = set(domain_df.loc[domain_df['pvalue'] < 0.05,'gene_id'].unique())
    intolerant_domains = set(domain_df.loc[domain_df['pvalue'] < 0.05, 'domain'].unique())
    
    # Split the DOMAINS column into multiple columns. Each domain is a column.
    domain_df = df["DOMAINS"].str.split('&', expand=True)
    # Identify all the ":" symbols in all domain column values and replace them with "_", using .replace() cannot
    # Replace ":" with "_" in all domain columns, handling both full and partial string matches
    # This will replace all instances of ":" even within substrings, e.g. "x:y:z" -> "x_y_z"
    for col in domain_df.columns:
        # Convert to string type first to handle any non-string values
        domain_df[col] = domain_df[col].astype(str)
        # Use str.replace with regex=False to replace all occurrences of ":" with "_"
        domain_df[col] = df["Gene"] + "_" + domain_df[col].str.replace(':', '_', regex=False)

    # Check if any domain in each row matches tolerant domains
    tolerant_domain_match = domain_df.isin(tolerant_domains).any(axis=1)
    intolerant_domain_match = domain_df.isin(intolerant_domains).any(axis=1)
    # Check if variant is missense and in a tolerant gene
    pp2_criteria = (df['Consequence'] == 'missense_variant') & (df['Gene'].isin(tolerant_genes) | tolerant_domain_match)
    # Check if variant is missense and in a intolerant domain
    bp1_criteria = (df['Consequence'] == 'missense_variant') & (df['Gene'].isin(intolerant_genes) | intolerant_domain_match)

    return pp2_criteria, bp1_criteria



def PP3_BP4_criteria(df: pd.DataFrame) -> pd.Series:
    # PP3: predicted to be deleterious by in-silico tools
    # Including PrimateAI (Missense), CADD, AlphaMissense (Missense), VEP, SpliceAI, SpliceVault, UTRAnnotator

    # BP4: variant is reported benign
    return (df['PrimateAI'] > 0.9) | \
           (df['CADD_phred'] >= 20) | \
           (df['am_class'] == 'pathogenic') | \
           (df['vep_consq_lof'] == True) | \
           (df['splicing_lof'] == True) | \
           (df['5UTR_lof'] == True), \
           (df['PrimateAI'] < 0.8) & \
           (df['CADD_phred'] < 20) & \
           (df['am_class'] != 'pathogenic') & \
           (df['vep_consq_lof'] != True) & \
           (df['splicing_lof'] != True) & \
           (df['5UTR_lof'] != True)


def PP5_BP6_criteria(df: pd.DataFrame) -> pd.Series:
    # PP5: The variant is reported as pathogenic by a reputable source but without to many supporting evidences
    high_confidence_status = {
        'practice_guideline',                                    # 4 stars
        'reviewed_by_expert_panel',                             # 3 stars
        'criteria_provided,_multiple_submitters,_no_conflicts'  # 2 stars
    }
    return df['CLNSIG'].str.contains('athogenic') & np.logical_not(df['CLNREVSTAT'].isin(high_confidence_status)), \
           df['CLNSIG'].str.contains('enign') & np.logical_not(df['CLNREVSTAT'].isin(high_confidence_status))



def BS1_criteria(df: pd.DataFrame, expected_incidence: float = 0.001) -> pd.Series:
    # BS1: PAF of variant is greater than expected incidence of the disease
    return df['gnomAD_joint_af'] > expected_incidence


def parse_hpo_inheritance(row_dict: dict) -> str:
    # Parse the HPO_gene_inheritance field and return the inheritance mode
    # The HPO_gene_inheritance field is a string with multiple inheritance modes separated by semicolons
    # These inheritance modes can correspond to 3 different pathogenic mechanisms: LoF, GoF, DN. 
    if isinstance(row_dict.get('HPO_gene_inheritance', None), str):
        hpo_inheritances = row_dict['HPO_gene_inheritance'].split(";")
    else:
        return None
    
    non_monogenic_set = {"Digenic inheritance", "Oligogenic inheritance", "Polygenic inheritance"}  # In most cases, these indicate compound heterozygous variants
    non_mendelian_set = {"Non-Mendelian inheritance"}  # Includes epigenetic modifications
    dominant_set = {"Autosomal dominant inheritance", "Autosomal dominant inheritance with maternal imprinting", "X-linked dominant inheritance"}
    recessive_set = {"Autosomal recessive inheritance", "X-linked recessive inheritance"}


    # HPO recessive
    hpo_recessive = all( hpo in recessive_set for hpo in hpo_inheritances )
    # HPO dominant
    hpo_dominant = all( hpo in dominant_set for hpo in hpo_inheritances )
    # HPO non monogenic
    hpo_non_monogenic = all( hpo in non_monogenic_set for hpo in hpo_inheritances )
    # HPO non mendelian
    hpo_non_mendelian = any( hpo in non_mendelian_set for hpo in hpo_inheritances )

    return {
            'hpo_recessive': hpo_recessive,
            'hpo_dominant': hpo_dominant,
            'hpo_non_monogenic': hpo_non_monogenic,
            'hpo_non_mendelian': hpo_non_mendelian
            }


def identify_inheritance_mode_per_row(row_dict: dict, gene_mean_am_score: float) -> Tuple[bool, bool]:
    # We need to use three fields of the table to determine the inheritance mode:
    # 1. Chromosome
    # 2. LOEUF
    # 3. AM score
    # 4. HPO_gene_inheritance (overrides the above two fields)

    haplo_insufficient = (float(row_dict.get('LOEUF', np.nan)) < 0.6) | (gene_mean_am_score > 0.58)
    haplo_sufficient = not haplo_insufficient

    hpo_inheritance = parse_hpo_inheritance(row_dict)
    if hpo_inheritance is None:
        return haplo_insufficient, haplo_sufficient

    dominant = False
    recessive = False

    if haplo_insufficient and \
       not any("recessive" in hpo for hpo in row_dict['HPO_gene_inheritance'].split(";")) and \
       not hpo_inheritance['hpo_non_mendelian'] and \
       not hpo_inheritance['hpo_non_monogenic']:
        # If haplo-insufficiency statistics indicate dominant inheritance, and HPO does not indicate recessive inheritance (nor non-mendelian and non-monogenic), then it is dominant
        dominant = True

    if haplo_sufficient and \
       not hpo_inheritance['hpo_dominant'] and \
       not hpo_inheritance['hpo_non_monogenic'] and \
       not hpo_inheritance['hpo_non_mendelian']:
        # If haplo-sufficiency statistics indicate recessive inheritance, and HPO does not indicate all inheritance found was dominant (nor non-mendelian and non-monogenic), then it is recessive
        # We can tolerate AD;AR simutaneously recorded in HPO for this gene because even for LoF it is haplo-sufficient, the AD can be for DN and GOF cases
        recessive = True

    if hpo_inheritance['hpo_recessive']:
        # If HPO indicates recessive inheritance only, then it is recessive
        recessive = True

    if hpo_inheritance['hpo_dominant']:
        # If HPO indicates dominant inheritance only, then it is dominant
        dominant = True

    return dominant, recessive

    

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
    dominant_array, recessive_array = zip(*results)
    return np.array(dominant_array), np.array(recessive_array)



def BS2_criteria(df: pd.DataFrame, 
                 gene_to_am_score_map: dict, 
                 threads: int = 10) -> pd.Series:
    '''
    The BS2 criteria is about observing variant in healthy adult
    There are several categories of situations here:
    1. For gene that is haplo-insufficient (autosomal dominant), we can assign BS2 if the variant is observed in a healthy adult (either homozygous or heterozygous)
    2. For gene that is autosomal recessive, we can assign BS2 if the variant is observed in a healthy adult with a homozygous genotype
    3. For X-linked recessive disease (gene on X but is haplo-sufficient), we can assign BS2 if the variant is observed in a healthy adult male (hemizygous) or a healthy adult female (homozygous)
    4. For X-linked dominant disease (gene on X and is haplo-insufficient), we can assign BS2 if the variant is observed in a healthy adult male (hemizygous)

    For haplo-insufficiency, we need to use LOEUF and mean AM score to determine.
    '''
    autosomal = (df['chrom'] != "chrX") & (df['chrom'] != "chrY")
    x_linked = df['chrom'] == "chrX"
    y_linked = df['chrom'] == "chrY"
    
    dominant, recessive = identify_inheritance_mode(df, gene_to_am_score_map, threads)

    # For autosomal dominant disease, we can assign BS2 if the variant is observed in a healthy adult (either homozygous or heterozygous)
    autosomal_dominant = autosomal & dominant & (df['gnomAD_joint_af'] > 0)
    autosomal_recessive = autosomal & recessive & ((df['gnomAD_nhomalt_XX'] > 0) | (df['gnomAD_nhomalt_XY'] > 0))

    # For X-linked disease, we can assign BS2 if the variant is observed in a healthy adult male (hemizygous) or a healthy adult female (homozygous)
    x_linked_recessive = x_linked & recessive & ((df['gnomAD_nhomalt_XX'] > 0) | (df['gnomAD_nhomalt_XY'] > 0))
    x_linked_dominant = x_linked & dominant & ((df['gnomAD_nhomalt_XX'] > 0) | (df['gnomAD_nhomalt_XY'] > 0) | (df['gnomAD_joint_af_XX'] > 0) | (df['gnomAD_joint_af_XY'] > 0))

    y_linked = y_linked & (df['gnomAD_nhomalt_XY'] > 0)

    return autosomal_dominant | autosomal_recessive | x_linked_recessive | x_linked_dominant | y_linked



def BS3_criteria(df: pd.DataFrame) -> pd.Series:
    # BS3: variant is reported benign and supported by functional evidences
    high_confidence_status = {
        'practice_guideline',                                    # 4 stars
        'reviewed_by_expert_panel',                             # 3 stars
        'criteria_provided,_multiple_submitters,_no_conflicts'  # 2 stars
    }
    return df['CLNSIG'].str.contains('enign') & df['CLNREVSTAT'].isin(high_confidence_status)



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


def BP2_PM3_criteria(df: pd.DataFrame, 
                     ped_df: pd.DataFrame, 
                     fam_name: str,
                     gene_to_am_score_map: dict,
                     ps1_criteria: pd.Series,
                     ps2_criteria: pd.Series,
                     ps3_criteria: pd.Series,
                     threads: int = 10) -> Tuple[pd.Series, pd.Series]:
    # BP2: observed in trans with a pathogenic variant in dominant disease, Or in-cis with a pathogenic variant with any inheritance mode
    # PM3: observed in trans with a pathogenic variant in recessive disease.
    pathogenic = df['vep_consq_lof'] | df['splicing_lof'] | df['5UTR_lof'] | ps1_criteria | (ps2_criteria & (df["CADD_phred"] >= 20)) | ps3_criteria
    is_dominant, is_recessive = identify_inheritance_mode(df, gene_to_am_score_map, threads)

    # Groupby gene and see if some variant is in-trans or in-cis with a pathogenic variant using the pathogenic boolean array above
    proband_info, father_info, mother_info, sib_info = identify_fam_members(ped_df, fam_name)
    proband, proband_pheno = proband_info
    
    all_genes = df['Gene'].unique().tolist()
    in_trans_pathogenic = np.array([False] * len(df))
    in_cis_pathogenic = np.array([False] * len(df))
    
    with Manager() as manager:
        shared_df = manager.dict()
        shared_df['data'] = df
        args = [(gene, shared_df['data'], pathogenic, proband) for gene in all_genes]

        with mp.Pool(threads) as pool:
            results = pool.starmap(check_gene_variants, args)
    
        for var_in_trans, var_in_cis in results:
            in_trans_pathogenic |= var_in_trans
            in_cis_pathogenic |= var_in_cis

    logger.info(f"There are {len(in_trans_pathogenic & is_dominant)} variants that are in-trans with pathogenic variants in a dominant (haploinsufficient) gene (disease)")
    logger.info(f"There are {len(in_trans_pathogenic & is_recessive)} variants that are in-trans with pathogenic variants in a recessive (haplo-sufficient) gene (disease)")
    logger.info(f"There are {len(in_cis_pathogenic)} variants that are in-cis with pathogenic variants regardless of the gene's known inheritance mode")

    bp2_criteria = (in_trans_pathogenic & is_dominant) | in_cis_pathogenic
    pm3_criteria = in_trans_pathogenic & is_recessive

    return bp2_criteria, pm3_criteria

    

def check_gene_variants(gene, df, pathogenic, proband):
    pathogenic_variants = df.loc[pathogenic, proband].tolist()
    var_in_trans = np.array([False] * len(df))
    var_in_cis = np.array([False] * len(df))
    if len([v for v in pathogenic_variants if v[-1] == "1"]) > 0:
        var_in_trans = ((df.loc[:, proband].str.split("|")[0] == "1") | (df.loc[:, proband].str.split("/")[0] == "1")) & (df.loc[:, "Gene"] == gene)
        var_in_cis = ((df.loc[:, proband].str.split("|")[1] == "1") | (df.loc[:, proband] == "1/1")) & (df.loc[:, "Gene"] == gene)
    if len([v for v in pathogenic_variants if v[0] == "1"]) > 0:
        var_in_trans = ((df.loc[:, proband].str.split("|")[1] == "1") | (df.loc[:, proband] == "1/1")) & (df.loc[:, "Gene"] == gene)
        var_in_cis = ((df.loc[:, proband].str.split("|")[0] == "1") | (df.loc[:, proband].str.split("/")[0] == "1")) & (df.loc[:, "Gene"] == gene)
    return var_in_trans, var_in_cis



def BP3_criteria(df: pd.DataFrame) -> pd.Series:
    # BP3: in-frame deletion in a repetitive region without a known function
    inframe_del = (df['Consequence'] == 'inframe_deletion') | \
                  (df['vep_consq_length_changing'] & ~df['vep_consq_lof']) | \
                  (df['splicing_len_changing'] & ~df['splicing_lof']) | \
                  (df['5UTR_len_changing'] & ~df['5UTR_lof'])
    repetitive_region = df['DOMAINS'].str.contains('Low_complexity')
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
        is_dominant, is_recessive = identify_inheritance_mode_per_row(row, gene_mean_am_score)
        
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
    no_splicing_altering = df['splicing_len_changing'] == False
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



def summary_acmg_per_var(row: pd.Series) -> pd.Series:
    '''
    Quantify ACMG criteria based on the ACMG_criteria column.
    
    Quantification:
    1. pvs = 9
    2. ps = 6
    3. pm = 2
    4. pp = 1.5
    
    1. bs = -6
    2. bp = -1.5
    
    Args:
        row: DataFrame row containing ACMG_criteria column
        
    Returns:
        Updated row with quantification results
    '''
    # Get criteria from ACMG_criteria column
    criteria = set(row['ACMG_criteria'].split(';')) if row['ACMG_criteria'] != '' else set()
    
    # Count criteria by type
    pvs_sum = sum(1 for c in criteria if c.startswith('PVS'))
    ps_sum = sum(1 for c in criteria if c.startswith('PS'))
    pm_sum = sum(1 for c in criteria if c.startswith('PM'))
    pp_sum = sum(1 for c in criteria if c.startswith('PP'))
    bs_sum = sum(1 for c in criteria if c.startswith('BS'))
    bp_sum = sum(1 for c in criteria if c.startswith('BP'))
    
    # Calculate final score
    final_score = pvs_sum * 9 + \
                  ps_sum * 6 + \
                  pm_sum * 2 + \
                  pp_sum * 1.5 + \
                  bs_sum * -6 + \
                  bp_sum * -1.5
    
    # Update row with results
    row['ACMG_quant_score'] = final_score
    
    # Assign final rank
    if final_score >= 12:
        row['ACMG_class'] = 'Pathogenic'
    elif final_score < 12 and final_score >= 6:
        row['ACMG_class'] = 'Likely_Pathogenic'
    elif final_score <= -12:
        row['ACMG_class'] = 'Benign'
    elif final_score <= -3 and final_score > -12:
        row['ACMG_class'] = 'Likely_Benign'
    else:
        row['ACMG_class'] = 'VOUS'
    
    return row



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



def ACMG_criteria_assign(anno_table: str, 
                         am_score_table: str, 
                         ped_table: str, 
                         fam_name: str,
                         clinvar_aa_dict_pkl: str,
                         intolerant_domains_pkl: str,
                         domain_mechanism_tsv: str,
                         alt_disease_vcf: str = "",
                         gnomAD_extreme_rare_threshold: float = 0.0001,
                         expected_incidence: float = 0.001,
                         threads: int = 10) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """
    Main function to assign ACMG criteria.
    Returns both annotated DataFrame and criteria matrix.
    """
    anno_df = pd.read_table(anno_table, low_memory=False)
    logger.info(f"Got {threads} threads to process the input table {anno_table}, now the table looks like: \n{anno_df[:5].to_string(index=False)}")
    anno_df = vep_consq_interpret(anno_df, threads)

    am_score_df = pd.read_table(am_score_table, low_memory=False)
    ped_df = pd.read_table(ped_table, low_memory=False)

    # Convert the am_score_df to a dictionary:
    # 1. Ensembl transcript ID (column 'transcript') to mean AM score (column 'mean_am_pathogenicity')
    am_score_dict = dict(zip(am_score_df['transcript'], am_score_df['mean_am_pathogenicity']))
    # Use anno_df to create a dict map from Ensembl transcript ID to gene ID
    transcript_to_gene_map = dict(zip(anno_df['Feature'], anno_df['Gene']))
    # Use the two dict above to create dict that maps gene ID to mean AM score
    gene_to_am_score_map = {g: am_score_dict[t] for t, g in transcript_to_gene_map.items() if t in am_score_dict}
    logger.info(f"gene_to_am_score_map created, {len(gene_to_am_score_map)} genes are having the AM score")

    # Apply the PVS1 criteria, LoF on a gene known to to be pathogenic due to LoF
    pvs1_criteria = PVS1_criteria(anno_df, gene_to_am_score_map)
    logger.info(f"PVS1 criteria applied, {pvs1_criteria.sum()} variants are having the PVS1 criteria")

    # Apply the PS1 and PM5 criteria
    ps1_criteria, pm5_criteria = PS1_PM5_criteria(anno_df, clinvar_aa_dict_pkl, threads)
    logger.info(f"PS1 criteria applied, {ps1_criteria.sum()} variants are having the PS1 criteria") 
    logger.info(f"PM5 criteria applied, {pm5_criteria.sum()} variants are having the PM5 criteria")

    # Apply the PS2 and PM6 criteria
    ps2_criteria, pm6_criteria = PS2_PM6_criteria(anno_df, ped_df, fam_name)
    logger.info(f"PS2 criteria applied, {ps2_criteria.sum()} variants are having the PS2 criteria")
    logger.info(f"PM6 criteria applied, {pm6_criteria.sum()} variants are having the PM6 criteria")

    # Apply the PS3 and BS3 criteria
    ps3bs3_results = PS3_BS3_criteria(anno_df)

    # Prevent double counting, if PVS1 is True, then PS3 should be False
    ps3_criteria = ps3bs3_results['PS3'] & ~pvs1_criteria
    bs3_criteria = ps3bs3_results['BS3']
    logger.info(f"PS3 criteria applied, {ps3_criteria.sum()} variants are having the PS3 criteria")
    logger.info(f"BS3 criteria applied, {bs3_criteria.sum()} variants are having the BS3 criteria")

    '''
    PS4 cannot be applied because usually we dont have enough cases to determine the frequency of the variant
    '''
    # Apply PM1 criteria, mutational hotspot or well-established functional protein domain
    pm1_criteria = PM1_criteria(anno_df, intolerant_domains_pkl, threads)
    logger.info(f"PM1 criteria applied, {pm1_criteria.sum()} variants are having the PM1 criteria")
    
    # Apply PM2 criteria, absent from gnomAD or extremely rare in gnomAD
    pm2_criteria = PM2_criteria(anno_df, gnomAD_extreme_rare_threshold)
    logger.info(f"PM2 criteria applied, {pm2_criteria.sum()} variants are having the PM2 criteria")
    
    # Apply PM4 criteria, causing the protein length change
    pm4_criteria = PM4_criteria(anno_df)
    
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

    # Apply PP3 criteria, predicted to be deleterious by in-silico tools
    pp3_criteria, bp4_criteria = PP3_BP4_criteria(anno_df)
    bp4_criteria = bp4_criteria & ~bs3_criteria
    logger.info(f"BP4 criteria applied, {bp4_criteria.sum()} variants are having the BP4 criteria")

    # Prevent double counting of PP3
    pp3_criteria = pp3_criteria & ~ps3_criteria & ~pvs1_criteria & ~pm4_criteria
    logger.info(f"PP3 criteria applied, {pp3_criteria.sum()} variants are having the PP3 criteria")
    '''
    PP4 cannot be applied, Patient
    's phenotype or family history is highly specific for a disease with a single genetic etiology
    '''
    # Apply PP5 criteria, reported as pathogenic by a reputable source but without to many supporting evidences
    # Apply BP6 criteria, reported as benign by a reputable source but without to many supporting evidences
    pp5_criteria, bp6_criteria = PP5_BP6_criteria(anno_df)

    '''
    BA1 cannot be applied, because usually not a single variant with PAF larger than 5% will survive to this step.
    '''
    # Apply BS1, PAF of variant is greater than expected incidence of the disease
    bs1_criteria = BS1_criteria(anno_df, expected_incidence)
    logger.info(f"BS1 criteria applied, {bs1_criteria.sum()} variants are having the BS1 criteria")

    # Apply BS2, variant observed in a healthy adult
    bs2_criteria = BS2_criteria(anno_df, gene_to_am_score_map, threads)
    logger.info(f"BS2 criteria applied, {bs2_criteria.sum()} variants are having the BS2 criteria")

    # Apply BS3, variant is reported benign and supported by functional evidences
    bs3_criteria = BS3_criteria(anno_df)
    logger.info(f"BS3 criteria applied, {bs3_criteria.sum()} variants are having the BS3 criteria")

    # Apply BS4, lack of family aggregation
    bs4_criteria = BS4_criteria(anno_df, ped_df, fam_name)
    logger.info(f"BS4 criteria applied, {bs4_criteria.sum()} variants are having the BS4 criteria")
    # Apply BP3, in-frame deletion in a repetitive region without a known function
    bp3_criteria = BP3_criteria(anno_df)
    logger.info(f"BP3 criteria applied, {bp3_criteria.sum()} variants are having the BP3 criteria")

    # Apply BP5, variant found in a sample with known alternative molecular basis for disease
    if alt_disease_vcf:
        bp5_criteria = BP5_criteria(anno_df, alt_disease_vcf, gene_to_am_score_map, threads)
    else:
        bp5_criteria = np.array([False] * len(anno_df))
    logger.info(f"BP5 criteria applied, {bp5_criteria.sum()} variants are having the BP5 criteria")

    # Apply BP7, synonymous variant, no splicing-altering consequence, not conserved. 
    bp7_criteria = BP7_criteria(anno_df)
    logger.info(f"BP7 criteria applied, {bp7_criteria.sum()} variants are having the BP7 criteria")
    # Apply BP2, observed in trans with a pathogenic variant in dominant disease, Or in-cis with a pathogenic variant with any inheritance mode
    # Apply PM3, observed in trans with a pathogenic variant in recessive disease.
    bp2_criteria, pm3_criteria = BP2_PM3_criteria(anno_df, 
                                                  ped_df, 
                                                  fam_name,
                                                  gene_to_am_score_map,
                                                  ps1_criteria,
                                                  ps2_criteria,
                                                  ps3_criteria,
                                                  threads)
    logger.info(f"BP2 criteria applied, {bp2_criteria.sum()} variants are having the BP2 criteria")
    logger.info(f"PM3 criteria applied, {pm3_criteria.sum()} variants are having the PM3 criteria")

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
    # Save the criteria matrix to a file which the path is based on the input anno_table
    output_matrix = ".".join(anno_table.split(".")[:-1]) + ".acmg.tsv"
    criteria_matrix.to_csv(output_matrix, sep="\t", index=False)
    
    # Apply quantification using ACMG_criteria column and use that to sort the variants
    anno_df = anno_df.apply(summary_acmg_per_var, axis=1)

    # Sort and rank variants
    anno_df = sort_and_rank_variants(anno_df)
    # Save the annotated table to replace the input anno_table
    anno_df.to_csv(anno_table, sep="\t", index=False)
    
    return anno_df, criteria_matrix



if __name__ == "__main__":
    parser = ap.ArgumentParser()
    parser.add_argument("--anno_table", type=str, required=True)
    parser.add_argument("--am_score_table", type=str, required=True)
    parser.add_argument("--ped_table", type=str, required=True)
    parser.add_argument("--fam_name", type=str, required=True)
    parser.add_argument("--clinvar_aa_dict_pkl", type=str, required=True)
    parser.add_argument("--intolerant_domains_pkl", type=str, required=True)
    parser.add_argument("--domain_mechanism_tsv", type=str, required=True)
    parser.add_argument("--alt_disease_vcf", type=str, required=False)
    parser.add_argument("--gnomAD_extreme_rare_threshold", type=float, required=False, default=0.0001)
    parser.add_argument("--expected_incidence", type=float, required=False, default=0.001)
    parser.add_argument("--threads", type=int, required=False, default=10)
    args = parser.parse_args()

    anno_df, criteria_matrix = ACMG_criteria_assign(args.anno_table, 
                                                    args.am_score_table, 
                                                    args.ped_table, 
                                                    args.fam_name,
                                                    args.clinvar_aa_dict_pkl,
                                                    args.intolerant_domains_pkl,
                                                    args.domain_mechanism_tsv,
                                                    alt_disease_vcf=args.alt_disease_vcf,
                                                    gnomAD_extreme_rare_threshold=args.gnomAD_extreme_rare_threshold,
                                                    expected_incidence=args.expected_incidence,
                                                    threads=args.threads)



