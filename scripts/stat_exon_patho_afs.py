#!/usr/bin/env python3

import os
import pysam
import numpy as np
import pickle
import argparse
import logging
import multiprocessing as mp
from collections import defaultdict
from functools import partial

# Set up logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger('vcf_exon_af_analyzer')

def extract_exon_from_csq(csq_data, transcript_field, exon_field):
    """Extract exon information from CSQ data."""
    exon_info = {}
    
    for csq in csq_data:
        fields = csq.split('|')
        if len(fields) <= max(transcript_field, exon_field):
            continue
            
        transcript = fields[transcript_field]
        exon = fields[exon_field]
            
        if transcript and exon:
            exon_info[transcript] = exon
            
    return exon_info

def process_chromosome(chromosome, vcf_path, csq_fields, high_confidence_status):
    """
    Process all variants on a chromosome and collect allele frequencies by exon.
    
    Args:
        chromosome: Chromosome name
        vcf_path: Path to VCF file
        csq_fields: Dictionary mapping field names to indices
        high_confidence_status: Set of high-confidence review statuses
        
    Returns:
        Dictionary mapping transcript to exon to list of allele frequencies
    """
    logger.info(f"Processing chromosome {chromosome}")
    
    # Create a nested defaultdict for exon-level allele frequencies
    transcript_exon_afs = {}
    
    try:
        vcf = pysam.VariantFile(vcf_path)
        count = 0
        
        # Fetch variants for this chromosome
        for record in vcf.fetch(chromosome):
            try:
                # Check if variant has clinical significance
                if 'CLNSIG' not in record.info or 'CLNREVSTAT' not in record.info or 'CSQ' not in record.info:
                    continue
                    
                # Check if pathogenic/likely pathogenic
                clnsig = ",".join(record.info['CLNSIG'])
                if not any(s in clnsig for s in ['Pathogenic', 'Likely_pathogenic']):
                    continue
                    
                # Check review status
                review_status = ",".join(record.info['CLNREVSTAT'])
                if not any(status in review_status for status in high_confidence_status):
                    continue
                
                # Get allele frequency
                af = None
                # Try different AF fields that might be in the VCF
                if "AF_grpmax_joint" in record.info:
                    af = record.info["AF_grpmax_joint"]
                    if isinstance(af, tuple):
                        af = float(af[0])  # Get first AF if it's a tuple

                if "AF_joint" in record.info and (af == 0 or af is None):
                    af = record.info["AF_joint"]
                    if isinstance(af, tuple):
                        af = float(af[0])  # Get first AF if it's a tuple
                
                if af is None or af <= 0:
                    af = 0
                
                # Parse transcript and exon information from CSQ
                csq_data = record.info['CSQ']
                exon_info = extract_exon_from_csq(
                    csq_data, 
                    csq_fields['Transcript'], 
                    csq_fields['EXON']
                )
                
                # Store allele frequency by transcript and exon
                for transcript, exon in exon_info.items():
                    if transcript and exon:
                        if transcript not in transcript_exon_afs:
                            transcript_exon_afs[transcript] = {}
                        if exon not in transcript_exon_afs[transcript]:
                            transcript_exon_afs[transcript][exon] = []
                        transcript_exon_afs[transcript][exon].append(af)
                        count += 1
                        
            except (KeyError, ValueError) as e:
                continue
        
        logger.info(f"Processed {count} pathogenic variants on chromosome {chromosome}")
        vcf.close()
        return transcript_exon_afs
        
    except Exception as e:
        logger.error(f"Error processing chromosome {chromosome}: {e}")
        return {}

def calculate_stats(transcript_exon_raw_afs):
    """
    Calculate statistics from raw allele frequencies.
    
    Args:
        transcript_exon_raw_afs: Dictionary mapping transcript to exon to list of AFs
        
    Returns:
        Dictionary with transcript -> exon -> (median, mean, max, min, np.array)
    """
    results = {}
    
    for transcript, exon_data in transcript_exon_raw_afs.items():
        results[transcript] = {}
        
        for exon, afs in exon_data.items():
            if not afs:
                continue
                
            afs_array = np.array(afs)
            
            # Calculate statistics
            median_af = float(np.median(afs_array))
            mean_af = float(np.mean(afs_array))
            max_af = float(np.max(afs_array))
            min_af = float(np.min(afs_array))
            
            # Store results
            results[transcript][exon] = (median_af, mean_af, max_af, min_af, afs_array)
    
    return results

def main():
    parser = argparse.ArgumentParser(
        description="Analyze allele frequencies of pathogenic variants by exon."
    )
    parser.add_argument('--vcf_path', required=True, help='Path to VCF file')
    parser.add_argument('--output', required=True, help='Output pickle file for results')
    parser.add_argument('--threads', type=int, required=False, default=4,
                        help='Number of processes to use (default: 4)')
    args = parser.parse_args()
    
    if args.threads is None:
        args.threads = 4
    
    # Define high-confidence statuses
    high_confidence_status = {
        'practice_guideline',                                    # 4 stars
        'reviewed_by_expert_panel',                              # 3 stars
        'criteria_provided,_multiple_submitters,_no_conflicts',  # 2 stars
    }
    
    try:
        # First pass to extract CSQ field structure and get chromosomes
        vcf = pysam.VariantFile(args.vcf_path)
        
        # Extract CSQ field indices
        csq_fields = {}
        if 'CSQ' in vcf.header.info:
            csq_header = vcf.header.info['CSQ'].description
            if 'Format:' in csq_header:
                csq_field_names = csq_header.split('Format: ')[1].split('|')
                for i, field_name in enumerate(csq_field_names):
                    csq_fields[field_name] = i
        
        if 'Transcript' not in csq_fields or 'EXON' not in csq_fields:
            # Try alternative field names
            if 'Feature' in csq_fields:
                csq_fields['Transcript'] = csq_fields['Feature']
            if 'EXON' not in csq_fields and 'exon' in csq_fields:
                csq_fields['EXON'] = csq_fields['exon']
                
        if 'Transcript' not in csq_fields or 'EXON' not in csq_fields:
            logger.error("Required CSQ fields not found in VCF header")
            return
        
        # Get list of chromosomes
        chromosomes = list(vcf.header.contigs)
        logger.info(f"Found {len(chromosomes)} chromosomes, processing in parallel with {args.threads} processes")
        vcf.close()
        
        # Process each chromosome in parallel
        with mp.Pool(args.threads) as pool:
            chromosome_results = pool.map(
                partial(
                    process_chromosome,
                    vcf_path=args.vcf_path,
                    csq_fields=csq_fields,
                    high_confidence_status=high_confidence_status
                ),
                chromosomes
            )
        
        # Combine results from all chromosomes
        all_transcript_exon_afs = defaultdict(lambda: defaultdict(list))
        for chrom_data in chromosome_results:
            for transcript, exon_data in chrom_data.items():
                for exon, afs in exon_data.items():
                    all_transcript_exon_afs[transcript][exon].extend(afs)
        
        # Calculate statistics for each transcript and exon
        final_results = calculate_stats(all_transcript_exon_afs)
        
        # Save results to pickle file
        with open(args.output, 'wb') as f:
            pickle.dump(final_results, f)
        
        # Report statistics
        total_transcripts = len(final_results)
        total_exons = sum(len(exons) for exons in final_results.values())
        logger.info(f"Analysis complete. Processed {total_transcripts} transcripts with {total_exons} exons.")
        logger.info(f"Results saved to {args.output}")
        
    except Exception as e:
        logger.error(f"Error: {e}")
        raise

if __name__ == "__main__":
    main()