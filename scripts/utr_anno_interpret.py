#!/usr/bin/env python
'''
Maintainer: yangyxt@hku.hk
This script is used to interpret the UTR annotations added by the UTRAnnotator (VEP plugin)
Basically there are these columns involved:
5UTR_annotation  5UTR_consequence  Existing_InFrame_oORFs  Existing_OutOfFrame_oORFs  Existing_uORFs
'''
import pandas as pd
import argparse as ap
from typing import Tuple, Dict, Set, Optional
import logging
import pickle
import multiprocessing as mp
from multiprocessing import shared_memory, Manager
import functools
import numpy as np
import time
import traceback
from io import StringIO
import sys
import uuid
import os

# Setup logging for main process
logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)
console_handler = logging.StreamHandler()
console_handler.setFormatter(logging.Formatter("%(levelname)s:%(asctime)s:%(funcName)s:%(lineno)s:%(message)s"))
logger.addHandler(console_handler)

# Log capture decorator for subprocess
def log_command(func):
    def wrapper(*args, **kwargs):
        # Create a string buffer to capture log output
        log_stream = StringIO()
        log_handler = logging.StreamHandler(log_stream)
        log_handler.setFormatter(logging.Formatter("%(levelname)s:%(asctime)s:%(funcName)s:%(lineno)s:%(message)s"))
        
        # Create subprocess logger
        row = args[0]  # First argument is row
        row_idx = row.name if hasattr(row, 'name') else 'unknown'
        subprocess_logger = logging.getLogger(f"subprocess-{row_idx}")
        subprocess_logger.setLevel(logging.INFO)
        subprocess_logger.addHandler(log_handler)
        
        # Call the function with the logger
        try:
            result = func(*args, logger=subprocess_logger, **kwargs)
            log_contents = log_stream.getvalue()
            return (True, result, log_contents)
        except Exception as e:
            tb_str = traceback.format_exc()
            subprocess_logger.error(f"Error: {e}\n{tb_str}")
            log_contents = log_stream.getvalue()
            return (False, (e, tb_str), log_contents)
    
    return wrapper

def interpret_uAUG_gained(annotation: str, logger=logger) -> Tuple[bool, bool]:
    """Interpret uAUG gained variants"""
    if not isinstance(annotation, str):
        return False, False
        
    fields = dict(item.split("=") for item in annotation.split(":"))
    confidence = (fields.get('KozakStrength') == 'Strong' or
                  (fields.get('KozakStrength') == 'Moderate' and 
                   fields.get('Evidence') == 'True'))
    
    is_lof = (
        fields.get('type') in ['OutOfFrame_oORF'] and
        confidence
    )
    
    # Length changing if it extends into CDS
    length_changing = (
        fields.get('type') in ['inFrame_oORF'] and
        confidence
    )

    logger.info(f"According to {annotation}, uAUG gained: {is_lof}, {length_changing}")
    
    return is_lof, length_changing or is_lof


def interpret_uSTOP_lost(annotation: str, logger=logger) -> Tuple[bool, bool]:
    """Interpret uSTOP lost variants"""
    if not isinstance(annotation, str):
        return False, False
        
    fields = dict(item.split("=") for item in annotation.split(":"))
    confidence = (fields.get('KozakStrength') == 'Strong' or
                  (fields.get('KozakStrength') == 'Moderate' and 
                   fields.get('Evidence') == 'True'))
    
    is_lof = (
        fields.get('AltStop') == 'False' and
        fields.get('FrameWithCDS') == 'outOfFrame' and
        confidence
    )
    
    length_changing = (
        confidence and 
        fields.get('AltStop') == 'False' and
        fields.get('FrameWithCDS') == 'inFrame'
    )
    
    logger.info(f"According to {annotation}, uSTOP lost: {is_lof}, {length_changing}")
    
    return is_lof, length_changing or is_lof


def interpret_uFrameshift(annotation: str, logger=logger) -> Tuple[bool, bool]:
    """Interpret uFrameshift variants"""
    if not isinstance(annotation, str):
        return False, False
        
    fields = dict(item.split("=") for item in annotation.split(":"))
    confidence = (fields.get('KozakStrength') == 'Strong' or
                  (fields.get('KozakStrength') == 'Moderate' and 
                   fields.get('Evidence') == 'True'))
    
    is_lof = (
        fields.get('alt_type_length') == 'NA' and
        fields.get('alt_type') == "OutOfFrame_oORF" and
        confidence
    )
    
    length_changing = (
        fields.get('alt_type') == "InFrame_oORF" and 
        fields.get('alt_type_length') == 'NA' and
        confidence
    )
    
    logger.info(f"According to {annotation}, uFrameshift: {is_lof}, {length_changing}")
    
    return is_lof, length_changing or is_lof


def check_domain_tolerance(transcript_id: str, 
                            exon_num: str,
                            domain_map: Dict,
                            intolerant_domains: Set[str],
                            logger=None) -> bool:
    """Check if length change affects intolerant domains"""
    if transcript_id not in domain_map:
        return False
        
    affected_domains = domain_map[transcript_id].get(exon_num, [])
    logger.info(f"Transcript {transcript_id} exon {exon_num} affected domains: {affected_domains}")

    if affected_domains:
        return any(domain in intolerant_domains for domain in affected_domains)
    else:
        return False


@log_command
def interpret_utr_annotation(row: pd.Series, 
                             domain_map: Optional[Dict] = None,
                             intolerant_domains: Optional[Set[str]] = None,
                             logger=None) -> Tuple[bool, bool]:
    """
    Main function to interpret UTR annotations
    
    Returns:
        Tuple[bool, bool]: (is_lof, length_changing)
    """
    row_idx = row.name
    logger.info(f"Processing row {row_idx}")
    consequence = row.get('5UTR_consequence', np.nan)
    annotation = row.get('5UTR_annotation', np.nan)
    
    if not isinstance(consequence, str) or not isinstance(annotation, str):
        logger.warning(f"Row {row_idx} has no consequence: {consequence} or annotation: {annotation}")
        return False, False
    
    # There can be multiple consequences and annotations for a variant, the values are separated by & symbol
    consequences = consequence.split('&')
    annotations = annotation.split('&')

    lofs, length_changes = [], []
    
    for consequence, annotation in zip(consequences, annotations):
        is_lof, length_changing = False, False
        # Updated consequence values based on UTRAnnotator documentation
        if '5_prime_UTR_premature_start_codon_gain_variant' in consequences:
            is_lof, length_changing = interpret_uAUG_gained(annotation, logger)
        elif '5_prime_UTR_uORF_stop_codon_loss_variant' in consequence:
            is_lof, length_changing = interpret_uSTOP_lost(annotation, logger)
        elif '5_prime_UTR_uORF_frameshift_variant' in consequence:
            is_lof, length_changing = interpret_uFrameshift(annotation, logger)
        
        # Check domain tolerance if length changing and domain info available
        if length_changing and domain_map and intolerant_domains:
            transcript_id = row.get('Feature', '')
            exon_num = "1"  # Since the variant is in the 5UTR, it must be in the first exon
            if check_domain_tolerance(transcript_id, exon_num, domain_map, intolerant_domains, logger):
                logger.info(f"Variant {row.get('chrom', '')}:{row.get('pos', '')}:{row.get('ref', '')}->{row.get('alt', '')} overlapping with {row.get('Feature', '')} affects exon 1, which overlaps with intolerant domains")
                is_lof = True
            
        lofs.append(is_lof)
        length_changes.append(length_changing)
        logger.info(f"Variant {row.get('chrom', '')}:{row.get('pos', '')}:{row.get('ref', '')}->{row.get('alt', '')} overlapping with {row.get('Feature', '')} has {consequence} consequence and {annotation} annotation. Is LOF: {is_lof}, Protein Length changing: {length_changing}")
    
    logger.info(f"Completed row {row_idx}: LOF={any(lofs)}, LengthChanging={any(length_changes)}")
    
    return any(lofs), any(length_changes)


# Define at module level
def process_row_shared(row_idx, shared_df_dict, shared_domain_map, shared_intolerant_domains):
    """Process a single row using shared data structures"""
    # Reconstruct the row from shared dict
    row_data = {col: shared_df_dict[col][row_idx] for col in shared_df_dict.keys()}
    row = pd.Series(row_data, name=row_idx)
    # Process the row with shared data
    return interpret_utr_annotation(row, shared_domain_map, shared_intolerant_domains)


def process_variants_table(variants_table: str,
                           domain_map_file: Optional[str] = None,
                           intolerant_domains_file: Optional[str] = None,
                           threads: Optional[int] = 12) -> pd.DataFrame:
    """Process entire variants table with parallel processing by rows"""
    
    df = pd.read_table(variants_table, low_memory=False)
    
    # Initialize result arrays
    lof_results = np.zeros(len(df), dtype=bool)
    len_changing_results = np.zeros(len(df), dtype=bool)
    
    # Create row indices
    row_indices = list(range(len(df)))
    
    # Load domain data in the main process
    domain_map = {}
    intolerant_domains = set()
    
    if domain_map_file and intolerant_domains_file:
        with open(domain_map_file, 'rb') as f:
            domain_map = pickle.load(f)
        with open(intolerant_domains_file, 'rb') as f:
            intolerant_domains = pickle.load(f)
    
    # Use a Manager to create shared data structures
    with Manager() as manager:
        # Create shared dictionary and set
        shared_domain_map = manager.dict(domain_map)
        shared_intolerant_domains = manager.list(intolerant_domains)
        
        # Also share the dataframe by converting it to a format the Manager can handle
        df_dict = {col: df[col].values for col in df.columns}
        shared_df_dict = manager.dict(df_dict)
        
        # Process rows in parallel using the shared data
        logger.info(f"Processing {len(df)} rows in parallel with {threads} threads")
        
        # Create a pool without initializer
        with mp.Pool(processes=threads) as pool:
            # Use partial to bind the shared data structures
            bound_process_row = functools.partial(
                process_row_shared,
                shared_df_dict=shared_df_dict,
                shared_domain_map=shared_domain_map,
                shared_intolerant_domains=shared_intolerant_domains
            )
            
            # Use map for ordered processing
            results = pool.imap(bound_process_row, row_indices)
            
            # Process results as before
            for i, (success, result, log_contents) in enumerate(results):
                print(f"\n************************************Row {i+1}/{len(df)} logs************************************\n", 
                     file=sys.stderr)
                print(log_contents, file=sys.stderr)
                print(f"\n************************************End Row {i+1}/{len(df)} logs************************************\n", 
                     file=sys.stderr)
                
                if success:
                    is_lof, len_changing = result
                    lof_results[i] = is_lof
                    len_changing_results[i] = len_changing
                else:
                    error_mes, tb_str = result
                    logger.error(f"Error processing row {i}: {error_mes}\n{tb_str}")
                    raise RuntimeError(f"Failed processing row {i}: {error_mes}")
    
    # Assign to dataframe
    df['5UTR_lof'] = lof_results
    df['5UTR_len_changing'] = len_changing_results
    
    df.to_csv(variants_table, sep='\t', index=False)
    return df


if __name__ == "__main__":
    parser = ap.ArgumentParser(description="Interpret UTR annotations")
    parser.add_argument("--variants_table", required=True, help="Path to the variants table")
    parser.add_argument("--domain_map", required=True, help="Path to the domain map, pickle file")
    parser.add_argument("--intolerant_domains", required=True, help="Path to the intolerant domains, pickle file")
    parser.add_argument("--threads", type=int, default=12, help="Number of processes to use")
    args = parser.parse_args()
    process_variants_table(args.variants_table, 
                           domain_map_file=args.domain_map, 
                           intolerant_domains_file=args.intolerant_domains,
                           threads=args.threads)

