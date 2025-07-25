#!/usr/bin/env python3
"""
Script to extract the most common pathogenic variant per gene from ClinVar VCF.
Uses pysam for parsing and multiprocessing for parallel processing.
Filters for high-confidence pathogenic variants based on review status.
Usage: python clinvar_most_common_pathogenic.py clinvar.vcf.gz output_pickle.pkl [num_threads]
"""

import sys
import pickle
import pysam
import re
from collections import defaultdict
import multiprocessing as mp
from functools import partial
import logging

# Configure logging
logging.basicConfig(level=logging.INFO, 
                    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

# Define high confidence review statuses
HIGH_CONFIDENCE_STATUS = {
    'practice_guideline',                                    # 4 stars
    'reviewed_by_expert_panel',                              # 3 stars
    'criteria_provided,_multiple_submitters,_no_conflicts',  # 2 stars
}

def get_chromosomes(vcf_file):
    """Get all chromosomes present in the VCF file"""
    with pysam.VariantFile(vcf_file) as vcf:
        return list(vcf.header.contigs)

def process_chromosome(chrom, vcf_file, csq_field_map):
    """Process variants on a single chromosome"""
    logger.info(f"Processing chromosome {chrom}")
    
    gene_to_variant = defaultdict(lambda: (None, 0))
    acc_patho_afs = defaultdict(lambda: 0)
    
    with pysam.VariantFile(vcf_file) as vcf:
        for record in vcf.fetch(chrom):
            # Skip if not in ClinVar
            if 'CLNSIG' not in record.info or 'CLNREVSTAT' not in record.info:
                continue
            
            # Get clinical significance
            clnsig = record.info.get('CLNSIG')
            if isinstance(clnsig, tuple):
                clnsig = ','.join(str(sig) for sig in clnsig)
            else:
                clnsig = str(clnsig)
                
            # Check if pathogenic or likely pathogenic
            if not ('athogenic' in clnsig):
                continue
            
            # Get review status
            clnrevstat = record.info.get('CLNREVSTAT')
            if isinstance(clnrevstat, tuple):
                clnrevstat = ','.join(str(stat) for stat in clnrevstat)
            else:
                clnrevstat = str(clnrevstat)
            
            # Check if high confidence review status
            is_high_confidence = False
            for status in HIGH_CONFIDENCE_STATUS:
                if status in clnrevstat.lower():
                    is_high_confidence = True
                    break
            
            if not is_high_confidence:
                continue
            
            # Extract gnomAD allele frequency
            af = record.info.get('AF_grpmax_joint', 0)
            if isinstance(af, tuple):
                af = af[0]

            if af is None or af == 0:
                af = record.info.get('AF_joint', 0)
                if isinstance(af, tuple):
                    af = af[0]

            if af is None or af <= 0:
                af = 0

            joint_af = record.info.get('AF_joint', 0)
            if isinstance(joint_af, tuple):
                joint_af = joint_af[0]

            if joint_af is None or joint_af <= 0:
                joint_af = 0
            
            # Skip if no CSQ field
            if 'CSQ' not in record.info:
                continue
            
            # Process each transcript annotation
            for csq in record.info['CSQ']:
                fields = csq.split('|')
                
                # Extract fields using the mapping
                try:
                    ensg = fields[csq_field_map['Gene']]
                    hgvsp = fields[csq_field_map['HGVSp']]
                    impact = fields[csq_field_map['IMPACT']]
                    consequence = fields[csq_field_map['Consequence']]
                    symbol = fields[csq_field_map['SYMBOL']]
                except (IndexError, KeyError):
                    continue
                
                if not ensg or not ensg.startswith('ENSG'):
                    continue
                
                # Check if this is a higher frequency than what we've seen before
                if ensg in acc_patho_afs:
                    acc_patho_afs[ensg] += joint_af
                else:
                    acc_patho_afs[ensg] = joint_af

                current_af = gene_to_variant[ensg][1]
                if af > current_af:
                    variant_info = {
                        'chrom': record.chrom,
                        'pos': record.pos,
                        'ref': record.ref,
                        'alt': ','.join(record.alts),
                        'af': af,
                        'hgvsp': hgvsp,
                        'impact': impact,
                        'consequence': consequence,
                        'clnsig': clnsig,
                        'clnrevstat': clnrevstat,
                        'gene_symbol': symbol,
                        'acc_patho_af': acc_patho_afs[ensg]
                    }
                    gene_to_variant[ensg] = (variant_info, af)
    
    # Convert to final dictionary format
    result = {gene: info for gene, (info, _) in gene_to_variant.items()}
    logger.info(f"Chromosome {chrom}: Found {len(result)} high-confidence pathogenic genes")
    return result

def get_csq_field_map(vcf_file):
    """Extract CSQ field format from VCF header and create a field name to index mapping"""
    with pysam.VariantFile(vcf_file) as vcf:
        for header_line in str(vcf.header).split('\n'):
            if 'ID=CSQ' in header_line:
                match = re.search(r'Format: ([\w\|\-]+)', header_line)
                if match:
                    fields = match.group(1).split('|')
                    return {field: idx for idx, field in enumerate(fields)}
    logger.error("Could not extract CSQ format from VCF header")
    return None

def merge_results(results):
    """Merge results from different chromosomes, keeping highest AF variant for each gene"""
    merged = {}
    for chrom_result in results:
        for gene, variant in chrom_result.items():
            if gene not in merged or merged[gene]['af'] < variant['af']:
                merged[gene] = variant
    return merged

def process_clinvar_vcf(vcf_file, output_pickle, num_threads=None):
    """Process ClinVar VCF and extract most common pathogenic variant per gene"""
    # Determine number of threads to use
    if num_threads is None:
        num_threads = min(mp.cpu_count(), 16)  # Default to CPU count but cap at 16
    else:
        num_threads = int(num_threads)
    
    logger.info(f"Processing {vcf_file} using {num_threads} threads")
    
    # Get CSQ field mapping
    csq_field_map = get_csq_field_map(vcf_file)
    if not csq_field_map:
        return
    
    # Get list of chromosomes
    chromosomes = get_chromosomes(vcf_file)
    logger.info(f"Found {len(chromosomes)} chromosomes: {', '.join(chromosomes[:5])}...")
    
    # Process each chromosome in parallel
    process_func = partial(process_chromosome, vcf_file=vcf_file, csq_field_map=csq_field_map)
    
    with mp.Pool(processes=num_threads) as pool:
        results = pool.map(process_func, chromosomes)
    
    # Merge results from all chromosomes
    final_results = merge_results(results)
    
    # Save to pickle
    with open(output_pickle, 'wb') as f:
        pickle.dump(final_results, f)
    
    logger.info(f"Processed {len(final_results)} genes with high-confidence pathogenic variants")
    logger.info(f"Results saved to {output_pickle}")
    
    return final_results

if __name__ == "__main__":
    if len(sys.argv) < 3:
        print(f"Usage: {sys.argv[0]} clinvar.vcf.gz output_pickle.pkl [num_threads]", file=sys.stderr)
        sys.exit(1)
    
    vcf_file = sys.argv[1]
    output_pickle = sys.argv[2]
    
    num_threads = None
    if len(sys.argv) > 3:
        num_threads = int(sys.argv[3])
    
    process_clinvar_vcf(vcf_file, output_pickle, num_threads)