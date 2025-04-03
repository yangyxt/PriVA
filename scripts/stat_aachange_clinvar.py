#!/usr/bin/env python3

import pysam
from collections import defaultdict
from typing import Dict, List, Tuple
import sys
import pickle
import logging
import multiprocessing as mp
from functools import partial

logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)
console_handler = logging.StreamHandler()
console_handler.setLevel(logging.INFO)
formatter = logging.Formatter("%(levelname)s:%(asctime)s:%(funcName)s:%(lineno)s:%(message)s")
console_handler.setFormatter(formatter)
logger.addHandler(console_handler)

# Define the nested defaultdict factory at the module level
def nested_defaultdict():
    return defaultdict(nested_defaultdict)

class AAClinvarCollector:
    def __init__(self, vcf_path: str, logger=logger):
        """Initialize collector with VCF file path."""
        self.vcf_path = vcf_path
        self.clinvaa_dict = nested_defaultdict()
        self.splice_dict = nested_defaultdict()  # New dict for splice-related data
        self.logger = logger
        self.high_confidence_status = {
            'practice_guideline',                                    # 4 stars
            'reviewed_by_expert_panel',                              # 3 stars
            'criteria_provided,_multiple_submitters,_no_conflicts',  # 2 stars
        }
        
    def get_csq_indices(self):
        """Extract CSQ format indices from VCF header."""
        vcf = pysam.VariantFile(self.vcf_path)
        csq_format = str(vcf.header).split('Format: ')[-1].strip('>"').split('|')
        
        # Get indices for all relevant fields
        indices = {
            'aa_pos_idx': csq_format.index('Protein_position') if 'Protein_position' in csq_format else -1,
            'hgvsp_idx': csq_format.index('HGVSp') if 'HGVSp' in csq_format else -1,
            'consq_idx': csq_format.index('Consequence') if 'Consequence' in csq_format else -1,
            'feature_type_idx': csq_format.index('Feature_type') if 'Feature_type' in csq_format else -1,
            'transcript_idx': csq_format.index('Feature') if 'Feature' in csq_format else -1,
            'exon_idx': csq_format.index('EXON') if 'EXON' in csq_format else -1,
            'intron_idx': csq_format.index('INTRON') if 'INTRON' in csq_format else -1,
            'hgvsc_idx': csq_format.index('HGVSc') if 'HGVSc' in csq_format else -1,
        }
        
        # Find SpliceAI field indices within CSQ format
        splice_ai_indices = {}
        for i, field in enumerate(csq_format):
            if field.startswith('SpliceAI_pred_'):
                splice_ai_indices[field] = i
                
        vcf.close()
        
        return indices, splice_ai_indices
    
    def get_chromosomes(self):
        """Get list of chromosomes from VCF file."""
        vcf = pysam.VariantFile(self.vcf_path)
        chromosomes = list(vcf.header.contigs)
        vcf.close()
        return chromosomes

    def is_pathogenic_with_good_review(self, clnsig, review_status):
        """Check if variant is pathogenic with good review status."""
        is_pathogenic = any(s in clnsig for s in ['Pathogenic', 'Likely_pathogenic'])
        has_good_review = any(status in review_status for status in self.high_confidence_status)
        return is_pathogenic and has_good_review

    def process_chromosome(self, chromosome: str, indices_and_fields: Tuple[Dict[str, int], Dict[str, int]]) -> Tuple[Dict, Dict]:
        """Process variants on a single chromosome, collecting both AA changes and splice data."""
        self.logger.info(f"Processing chromosome {chromosome}")
        
        indices, splice_ai_indices = indices_and_fields
        
        # Use regular dicts instead of nested_defaultdict for pickling
        aa_chrom_results = {}
        splice_chrom_results = {}
        
        try:
            vcf = pysam.VariantFile(self.vcf_path)
            count_aa = 0
            count_splice = 0
            
            # Extract indices
            aa_pos_idx = indices['aa_pos_idx']
            hgvsp_idx = indices['hgvsp_idx']
            consq_idx = indices['consq_idx']
            feature_type_idx = indices['feature_type_idx']
            transcript_idx = indices['transcript_idx']
            exon_idx = indices['exon_idx']
            intron_idx = indices['intron_idx']
            hgvsc_idx = indices['hgvsc_idx']
            
            # Process each variant on this chromosome
            for record in vcf.fetch(chromosome):
                try:
                    if 'CLNREVSTAT' not in record.info or 'CLNSIG' not in record.info:
                        continue

                    rev_status = ",".join(record.info['CLNREVSTAT'])
                    cln_sig = ",".join(record.info['CLNSIG'])
                    
                    # Only process pathogenic/likely pathogenic variants with good review status
                    if not self.is_pathogenic_with_good_review(cln_sig, rev_status):
                        continue

                    csq_value = record.info['CSQ']
                    if isinstance(csq_value, tuple):
                        transcript_annotations = csq_value
                    else:
                        transcript_annotations = csq_value.split(',')

                    for annotation in transcript_annotations:
                        fields = annotation.split('|')

                        # Skip if fields are incomplete
                        if len(fields) <= max(
                            feature_type_idx, transcript_idx, 
                            aa_pos_idx, hgvsp_idx, consq_idx, 
                            exon_idx, intron_idx, hgvsc_idx
                        ):
                            continue

                        enst_id = fields[transcript_idx]
                        if not enst_id:
                            continue
                        
                        # Process amino acid changes
                        aa_pos = fields[aa_pos_idx] if aa_pos_idx >= 0 else None
                        hgvsp = fields[hgvsp_idx] if hgvsp_idx >= 0 else None
                        
                        if aa_pos and hgvsp:
                            # Create nested structure if needed for AA changes
                            if enst_id not in aa_chrom_results:
                                aa_chrom_results[enst_id] = {}
                            if aa_pos not in aa_chrom_results[enst_id]:
                                aa_chrom_results[enst_id][aa_pos] = {}
                            if hgvsp not in aa_chrom_results[enst_id][aa_pos]:
                                aa_chrom_results[enst_id][aa_pos][hgvsp] = {}
                                
                            current_dict = aa_chrom_results[enst_id][aa_pos][hgvsp]

                            if "CLNSIG" not in current_dict:
                                current_dict["CLNSIG"] = [cln_sig]
                            else:
                                current_dict["CLNSIG"].append(cln_sig)

                            if "CLNREVSTAT" not in current_dict:
                                current_dict["CLNREVSTAT"] = [rev_status]
                            else:
                                current_dict["CLNREVSTAT"].append(rev_status)
                                
                            count_aa += 1
                        
                        # Process splice-related data
                        exon = fields[exon_idx] if exon_idx >= 0 else None
                        intron = fields[intron_idx] if intron_idx >= 0 else None
                        hgvsc = fields[hgvsc_idx] if hgvsc_idx >= 0 else None
                        consq = fields[consq_idx] if consq_idx >= 0 else None
                        
                        # Extract SpliceAI values from CSQ fields
                        splice_ai_data = {}
                        for field_name, idx in splice_ai_indices.items():
                            if idx < len(fields) and fields[idx]:
                                splice_ai_data[field_name] = fields[idx]
                        
                        # Continue only if we have at least one piece of splice-related data
                        # logger.info(f"The record has exon {exon}, intron {intron}, consq {consq}, rev_status {rev_status}, and cln_sig {cln_sig}")
                        if (exon or intron) and ("athogenic" in cln_sig) and (rev_status in self.high_confidence_status):
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
                            
                            # Create nested structure for splice data
                            if enst_id not in splice_chrom_results:
                                splice_chrom_results[enst_id] = []
                            
                            splice_chrom_results[enst_id].append(pos_info)
                            count_splice += 1

                except (KeyError, ValueError) as e:
                    self.logger.warning(f"Error processing record {record.chrom}:{record.pos} - {str(e)}")
                    continue
                    
            self.logger.info(f"Processed {count_aa} AA changes and {count_splice} splice variants on chromosome {chromosome}")
            vcf.close()
            return (aa_chrom_results, splice_chrom_results)
            
        except Exception as e:
            self.logger.error(f"Error processing chromosome {chromosome}: {e}")
            return ({}, {})

    def collect_data(self, n_processes=None):
        """Process VCF file and collect both AA changes and splice data in parallel."""
        if n_processes is None:
            n_processes = 4
            
        self.logger.info(f"Processing with {n_processes} processes")
        
        # Get CSQ indices and SpliceAI fields once
        indices_and_fields = self.get_csq_indices()
        
        # Get list of chromosomes
        chromosomes = self.get_chromosomes()
        self.logger.info(f"Found {len(chromosomes)} chromosomes")
        
        # Process chromosomes in parallel
        with mp.Pool(n_processes) as pool:
            chrom_results = pool.map(
                partial(self.process_chromosome, indices_and_fields=indices_and_fields),
                chromosomes
            )
            
        # Combine results into the nested defaultdict structures
        aa_total = 0
        splice_total = 0
        for aa_data, splice_data in chrom_results:
            aa_count = self.merge_aa_data(aa_data)
            splice_count = self.merge_splice_data(splice_data)
            aa_total += aa_count
            splice_total += splice_count
            
        self.logger.info(f"Finished collecting data: {aa_total} AA changes and {splice_total} splice variants")

    def merge_aa_data(self, aa_data):
        """Merge AA data into the main nested defaultdict."""
        count = 0
        for enst_id, aa_pos_data in aa_data.items():
            for aa_pos, hgvsp_data in aa_pos_data.items():
                for hgvsp, info in hgvsp_data.items():
                    current_dict = self.clinvaa_dict[enst_id][aa_pos][hgvsp]
                    
                    if "CLNSIG" not in current_dict and "CLNSIG" in info:
                        current_dict["CLNSIG"] = info["CLNSIG"]
                    elif "CLNSIG" in info:
                        current_dict["CLNSIG"].extend(info["CLNSIG"])
                        
                    if "CLNREVSTAT" not in current_dict and "CLNREVSTAT" in info:
                        current_dict["CLNREVSTAT"] = info["CLNREVSTAT"]
                    elif "CLNREVSTAT" in info:
                        current_dict["CLNREVSTAT"].extend(info["CLNREVSTAT"])
                    count += 1
        return count

    def merge_splice_data(self, splice_data):
        """Merge splice data into the splice nested defaultdict."""
        count = 0
        for enst_id, variants in splice_data.items():
            if enst_id not in self.splice_dict:
                self.splice_dict[enst_id] = []
            self.splice_dict[enst_id].extend(variants)
            count += len(variants)
        return count

    # Keep original method for backward compatibility
    def collect_clinvar_aa(self, n_processes=None):
        """Legacy method that calls collect_data."""
        self.collect_data(n_processes)

    @staticmethod
    def load_scores_from_pickle(pickle_path: str) -> Dict:
        """Load scores from the pickle file."""
        with open(pickle_path, 'rb') as f:
            data = pickle.load(f)
        return data

def main(vcf_path, output_pickle, splice_output_pickle=None, n_processes=None):
    """
    Notice that the input VCF file must be the ClinVar VCF file annotated by VEP.
    
    Args:
        vcf_path: Path to VCF file
        output_pickle: Path to output pickle file for AA changes
        splice_output_pickle: Path to output pickle file for splice data (optional)
        n_processes: Number of processes to use (default: CPU count)
    """
    if n_processes is None:
        n_processes = 4 # default number of processes
        
    collector = AAClinvarCollector(vcf_path)
    collector.collect_data(n_processes=n_processes)

    # Output the AA changes to a pickle file
    if output_pickle:
        with open(output_pickle, 'wb') as f:
            pickle.dump(collector.clinvaa_dict, f)
        logger.info(f"AA change results saved to {output_pickle}")
    
    # Output the splice data to a pickle file if specified
    if splice_output_pickle:
        with open(splice_output_pickle, 'wb') as f:
            pickle.dump(collector.splice_dict, f)
        logger.info(f"Splice data results saved to {splice_output_pickle}")

    return collector.clinvaa_dict, collector.splice_dict

if __name__ == "__main__":
    if len(sys.argv) < 3:
        print("Usage: python stat_aachange_clinvar.py <vcf_path> <aa_output_pickle> [splice_output_pickle] [num_processes]")
        sys.exit(1)
        
    vcf_path = sys.argv[1]
    output_pickle = sys.argv[2]
    
    splice_output_pickle = None
    if len(sys.argv) > 3:
        splice_output_pickle = sys.argv[3]
    
    n_processes = None
    if len(sys.argv) > 4:
        n_processes = int(sys.argv[4])
        
    main(vcf_path, output_pickle, splice_output_pickle, n_processes)