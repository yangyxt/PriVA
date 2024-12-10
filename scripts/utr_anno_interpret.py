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

logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)
console_handler = logging.StreamHandler()
console_handler.setFormatter(logging.Formatter("%(levelname)s:%(asctime)s:%(funcName)s:%(lineno)s:%(message)s"))
logger.addHandler(console_handler)


def interpret_uAUG_gained(annotation: str) -> Tuple[bool, bool]:
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
    
    return is_lof, length_changing or is_lof


def interpret_uSTOP_lost(annotation: str) -> Tuple[bool, bool]:
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
    
    return is_lof, length_changing or is_lof



def interpret_uFrameshift(annotation: str) -> Tuple[bool, bool]:
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
    
    return is_lof, length_changing or is_lof



def check_domain_tolerance(transcript_id: str, 
                         exon_num: str,
                         domain_map: Dict,
                         intolerant_domains: Set[str]) -> bool:
    """Check if length change affects intolerant domains"""
    if transcript_id not in domain_map:
        return False
        
    affected_domains = domain_map[transcript_id].get(exon_num, [])
    logger.info(f"Transcript {transcript_id} exon {exon_num} affected domains: {affected_domains}")
    return any(domain in intolerant_domains for domain in affected_domains)



def interpret_utr_annotation(row: pd.Series, 
                           domain_map: Optional[Dict] = None,
                           intolerant_domains: Optional[Set[str]] = None) -> Tuple[bool, bool]:
    """
    Main function to interpret UTR annotations
    
    Returns:
        Tuple[bool, bool]: (is_lof, length_changing)
    """
    consequence = row.get('5UTR_consequence', '')
    annotation = row.get('5UTR_annotation', '')
    
    if not isinstance(consequence, str) or not isinstance(annotation, str):
        return False, False
    
    # There can be multiple consequences and annotations for a variant, the values are separated by & symbol
    consequences = consequence.split('&')
    annotations = annotation.split('&')

    lofs, length_changes = [], []
    
    for consequence, annotation in zip(consequences, annotations):
        is_lof, length_changing = False, False
        # Updated consequence values based on UTRAnnotator documentation
        if '5_prime_UTR_premature_start_codon_gain_variant' in consequences:
            is_lof, length_changing = interpret_uAUG_gained(annotation)
        elif '5_prime_UTR_uORF_stop_codon_loss_variant' in consequence:
            is_lof, length_changing = interpret_uSTOP_lost(annotation)
        elif '5_prime_UTR_uORF_frameshift_variant' in consequence:
            is_lof, length_changing = interpret_uFrameshift(annotation)
        
        # Check domain tolerance if length changing and domain info available
        if length_changing and domain_map and intolerant_domains:
            transcript_id = row.get('Feature', '')
            exon_num = "1"  # Since the variant is in the 5UTR, it must be in the first exon
            if check_domain_tolerance(transcript_id, exon_num, domain_map, intolerant_domains):
                logger.info(f"Variant {row.chrom}:{row.pos}:{row.ref}->{row.alt} overlapping with {row.Feature} affects exon 1, which overlaps with intolerant domains")
                is_lof = True
            
        lofs.append(is_lof)
        length_changes.append(length_changing)
        logger.info(f"Variant {row.chrom}:{row.pos}:{row.ref}->{row.alt} overlapping with {row.Feature} has {consequence} consequence and {annotation} annotation. Is LOF: {is_lof}, Protein Length changing: {length_changing}")
    
    return any(lofs), any(length_changes)



def process_row(args) -> Tuple[bool, bool]:
    """Process a single row with domain information"""
    row, domain_map, intolerant_domains = args
    is_lof, length_changing = interpret_utr_annotation(
        row, domain_map, intolerant_domains
    )
    return is_lof, length_changing



def process_variants_table(variants_table: str,
                         domain_map_file: Optional[str] = None,
                         intolerant_domains_file: Optional[str] = None,
                         threads: Optional[int] = 12) -> pd.DataFrame:
    """Process entire variants table with parallel processing by rows"""
    
    df = pd.read_table(variants_table, low_memory=False)
    
    if domain_map_file and intolerant_domains_file:
        domain_map = pickle.load(open(domain_map_file, 'rb'))
        intolerant_domains = pickle.load(open(intolerant_domains_file, 'rb'))
    
    # Create arguments for each row
    row_args = [(row, domain_map, intolerant_domains) for _, row in df.iterrows()]
    
    # Process rows in parallel
    with mp.Pool(processes=threads) as pool:
        results = pool.map(process_row, row_args)
    
    df['5UTR_lof'], df['5UTR_len_changing'] = zip(*results)
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

