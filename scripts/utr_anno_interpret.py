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
import numpy as np
import time
import sys

# Setup logging for main process
logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)
console_handler = logging.StreamHandler()
console_handler.setFormatter(logging.Formatter("%(levelname)s:%(asctime)s:%(funcName)s:%(lineno)s:%(message)s"))
logger.addHandler(console_handler)

def interpret_uAUG_gained(annotation: str, logger=logger) -> Tuple[bool, bool]:
    """Interpret uAUG gained variants"""
    if not isinstance(annotation, str):
        return False, False, False
        
    fields = dict(item.split("=") for item in annotation.split(":"))
    confidence = (fields.get('KozakStrength') == 'Strong' or
                  (fields.get('KozakStrength') == 'Moderate' and 
                   fields.get('Evidence') == 'True'))
    
    offset = int(fields.get('DistanceToCDS', 0))
    
    is_lof = (
        fields.get('type') in ['OutOfFrame_oORF'] and
        confidence
    )

    frameshift = (
        fields.get('type') in ['OutOfFrame_oORF'] and
        confidence
    )
    
    # Length changing if it extends into CDS
    length_changing = (
        fields.get('type') in ['inFrame_oORF'] and
        confidence
    )

    logger.info(f"According to {annotation}, uAUG gained: {is_lof}, {length_changing}")
    
    return is_lof, frameshift, length_changing or is_lof


def interpret_uSTOP_lost(annotation: str, logger=logger) -> Tuple[bool, bool]:
    """Interpret uSTOP lost variants"""
    if not isinstance(annotation, str):
        return False, False, False
        
    fields = dict(item.split("=") for item in annotation.split(":"))
    confidence = (fields.get('KozakStrength') == 'Strong' or
                  (fields.get('KozakStrength') == 'Moderate' and 
                   fields.get('Evidence') == 'True'))
    
    is_lof = (
        fields.get('AltStop') == 'False' and
        fields.get('FrameWithCDS') == 'outOfFrame' and
        confidence
    )

    frameshift = (
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
    
    return is_lof, frameshift, length_changing or is_lof


def interpret_uFrameshift(annotation: str, logger=logger) -> Tuple[bool, bool]:
    """Interpret uFrameshift variants"""
    if not isinstance(annotation, str):
        return False, False, False
        
    fields = dict(item.split("=") for item in annotation.split(":"))
    confidence = (fields.get('KozakStrength') == 'Strong' or
                  (fields.get('KozakStrength') == 'Moderate' and 
                   fields.get('Evidence') == 'True'))
    
    is_lof = (
        fields.get('alt_type_length') == 'NA' and
        fields.get('alt_type') == "OutOfFrame_oORF" and
        confidence
    )

    frameshift = (
        fields.get('alt_type') == "OutOfFrame_oORF" and
        fields.get('alt_type_length') == 'NA' and
        confidence
    )
    
    length_changing = (
        fields.get('alt_type') == "InFrame_oORF" and 
        fields.get('alt_type_length') == 'NA' and
        confidence
    )
    
    logger.info(f"According to {annotation}, uFrameshift: {is_lof}, {length_changing}")
    
    return is_lof, frameshift, length_changing or is_lof


def check_domain_tolerance(transcript_id: str, 
                           exon_num: str,
                           domain_map: Dict,
                           intolerant_domains: Set[str],
                           logger=logger) -> bool:
    """Check if length change affects intolerant domains"""
    if transcript_id not in domain_map:
        return False
        
    affected_domains = domain_map[transcript_id].get(exon_num, [])
    logger.info(f"Transcript {transcript_id} exon {exon_num} affected domains: {affected_domains}")

    if affected_domains:
        return any(domain in intolerant_domains for domain in affected_domains)
    else:
        return False


def interpret_utr_annotation(row: pd.Series, 
                            domain_map: Optional[Dict] = None,
                            intolerant_domains: Optional[Set[str]] = None) -> Tuple[bool, bool]:
    """
    Main function to interpret UTR annotations
    
    Returns:
        Tuple[bool, bool]: (is_lof, length_changing)
    """
    row_idx = row.name
    logger.debug(f"Processing row {row_idx}")
    consequence = row.get('5UTR_consequence', '')
    annotation = row.get('5UTR_annotation', '')
    
    if not isinstance(consequence, str) or not isinstance(annotation, str):
        logger.debug(f"Row {row_idx} has no consequence: {consequence} or annotation: {annotation}")
        return False, False, False, False
    
    # There can be multiple consequences and annotations for a variant, the values are separated by & symbol
    consequences = consequence.split('&')
    annotations = annotation.split('&')

    lofs, frameshifts, length_changes, span_intol_domains = [], [], [], []
    
    for consequence, annotation in zip(consequences, annotations):
        is_lof, frameshift, length_changing, span_intol_domain = False, False, False, False
        # Updated consequence values based on UTRAnnotator documentation
        if '5_prime_UTR_premature_start_codon_gain_variant' in consequences:
            is_lof, frameshift, length_changing = interpret_uAUG_gained(annotation)
        elif '5_prime_UTR_uORF_stop_codon_loss_variant' in consequence:
            is_lof, frameshift, length_changing = interpret_uSTOP_lost(annotation)
        elif '5_prime_UTR_uORF_frameshift_variant' in consequence:
            is_lof, frameshift, length_changing = interpret_uFrameshift(annotation)
        
        # Check domain tolerance if length changing and domain info available
        if length_changing and domain_map and intolerant_domains:
            transcript_id = row.get('Feature', '')
            exon_num = "1"  # Since the variant is in the 5UTR, it must be in the first exon
            if check_domain_tolerance(transcript_id, exon_num, domain_map, intolerant_domains):
                logger.debug(f"Variant {row.get('chrom', '')}:{row.get('pos', '')}:{row.get('ref', '')}->{row.get('alt', '')} overlapping with {row.get('Feature', '')} affects exon 1, which overlaps with intolerant domains")
                is_lof = True
                span_intol_domain = True
            
        lofs.append(is_lof)
        frameshifts.append(frameshift)
        length_changes.append(length_changing)
        span_intol_domains.append(span_intol_domain)
        logger.debug(f"Variant {row.get('chrom', '')}:{row.get('pos', '')}:{row.get('ref', '')}->{row.get('alt', '')} overlapping with {row.get('Feature', '')} has {consequence} consequence and {annotation} annotation. Is LOF: {is_lof}, Protein Length changing: {length_changing}")
    
    logger.debug(f"Completed row {row_idx}: LOF={any(lofs)}, LengthChanging={any(length_changes)}")
    
    return any(lofs), any(frameshifts), any(length_changes), any(span_intol_domains)


def process_variants_table(variants_table: str,
                           domain_map_file: Optional[str] = None,
                           intolerant_domains_file: Optional[str] = None) -> pd.DataFrame:
    """Process variants table with serial processing"""
    
    start_time = time.time()
    logger.info(f"Loading variants table from {variants_table}")
    df = pd.read_table(variants_table, low_memory=False)
    logger.info(f"Loaded {len(df)} variants in {time.time() - start_time:.2f} seconds")
    
    # Initialize result arrays
    lof_results = np.zeros(len(df), dtype=bool)
    frameshift_results = np.zeros(len(df), dtype=bool)
    len_changing_results = np.zeros(len(df), dtype=bool)
    span_intol_domain_results = np.zeros(len(df), dtype=bool)
    
    # Load domain data once
    domain_map = {}
    intolerant_domains = set()
    
    if domain_map_file and intolerant_domains_file:
        logger.info(f"Loading domain data from {domain_map_file} and {intolerant_domains_file}")
        start_time = time.time()
        with open(domain_map_file, 'rb') as f:
            domain_map = pickle.load(f)
        with open(intolerant_domains_file, 'rb') as f:
            intolerant_domains = pickle.load(f)
        logger.info(f"Loaded domain data in {time.time() - start_time:.2f} seconds")
    
    # Process each row serially
    logger.info(f"Processing {len(df)} rows serially")
    start_time = time.time()
    processed = 0
    report_interval = max(1, min(10000, len(df) // 20))  # Report progress every ~5%
    
    for idx in range(len(df)):
        row = df.iloc[idx]
        is_lof, frameshift, len_changing, span_intol_domain = interpret_utr_annotation(row, domain_map, intolerant_domains)
        lof_results[idx] = is_lof
        frameshift_results[idx] = frameshift
        len_changing_results[idx] = len_changing
        span_intol_domain_results[idx] = span_intol_domain
        
        processed += 1
        
        # Report progress periodically
        if processed % report_interval == 0:
            elapsed = time.time() - start_time
            rate = processed / elapsed if elapsed > 0 else 0
            estimated_total = len(df) / rate if rate > 0 else 0
            remaining = estimated_total - elapsed
            
            logger.info(f"Processed {processed}/{len(df)} rows ({processed/len(df)*100:.1f}%) "
                      f"at {rate:.1f} rows/sec. Estimated time remaining: {remaining/60:.1f} minutes")
    
    total_time = time.time() - start_time
    logger.info(f"Completed processing {len(df)} rows in {total_time:.2f} seconds "
               f"({len(df)/total_time:.1f} rows/sec)")
    
    # Assign to dataframe
    df['5UTR_lof'] = lof_results
    df['5UTR_len_changing'] = len_changing_results
    df['5UTR_frameshift'] = frameshift_results
    df['5UTR_span_intol_domain'] = span_intol_domain_results
    
    logger.info(f"Saving results to {variants_table}")
    df.to_csv(variants_table, sep='\t', index=False)
    return df


if __name__ == "__main__":
    parser = ap.ArgumentParser(description="Interpret UTR annotations")
    parser.add_argument("--variants_table", required=True, help="Path to the variants table")
    parser.add_argument("--domain_map", required=True, help="Path to the domain map, pickle file")
    parser.add_argument("--intolerant_domains", required=True, help="Path to the intolerant domains, pickle file")
    args = parser.parse_args()
    
    process_variants_table(args.variants_table, 
                         domain_map_file=args.domain_map, 
                         intolerant_domains_file=args.intolerant_domains)

