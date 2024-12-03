#!/usr/bin/env python3

import pysam
import numpy as np
from collections import defaultdict
from typing import Dict, List
import sys
import pickle

# Define the nested defaultdict factory at the module level
def nested_defaultdict():
    return defaultdict(nested_defaultdict)

class DomainAMScoreCollector:
    def __init__(self, vcf_path: str):
        """Initialize collector with VCF file path."""
        self.vcf_path = vcf_path
        self.scores_dict = nested_defaultdict()
        
    def _parse_domain_hierarchy(self, domains_str: str) -> List[tuple]:
        """Parse domain string into hierarchical components.
        
        Args:
            domains_str: String like 'PANTHER:PTHR26451&PANTHER:PTHR26451:SF72'
            
        Returns:
            List of tuples containing (db_name, hierarchy_list)
            e.g., [('PANTHER', ['PTHR26451']), ('PANTHER', ['PTHR26451', 'SF72'])]
        """
        if not domains_str or domains_str == '':
            return []
            
        domain_entries = domains_str.split('&')
        parsed_domains = []
        
        for entry in domain_entries:
            parts = entry.split(':')
            db_name = parts[0]
            
            # Build hierarchical structure
            for i in range(1, len(parts)):
                hierarchy = parts[1:i+1]  # Include all levels up to current
                parsed_domains.append((db_name, hierarchy))
                
        return parsed_domains

    def collect_scores(self):
        """Process VCF file and collect AlphaMissense pathogenicity scores."""
        vcf = pysam.VariantFile(self.vcf_path)
        
        # Get CSQ format from header
        csq_format = str(vcf.header).split('Format: ')[-1].strip('>"').split('|')
        domains_idx = csq_format.index('DOMAINS')
        gene_idx = csq_format.index('Gene')
        exon_idx = csq_format.index('EXON')
        feature_type_idx = csq_format.index('Feature_type')
        transcript_idx = csq_format.index('Feature')
        
        for record in vcf:
            try:
                am_score = float(record.info['AM_PATHOGENICITY'])
                prot_var = record.info['PVAR']
                # Extract amino acid position from PVAR (assuming format like 'p.Arg123Ser')
                aa_pos = prot_var[:-1]

                csq_value = record.info['CSQ']
                if isinstance(csq_value, tuple):
                    transcript_annotations = csq_value
                else:
                    transcript_annotations = csq_value.split(',')

                for csq_entry in transcript_annotations:
                    fields = csq_entry.split('|')
                    
                    # Skip if not a transcript or not an ENST transcript
                    if (fields[feature_type_idx] != 'Transcript' or 
                        not fields[transcript_idx].startswith('ENST')):
                        continue
                    
                    domains_str = fields[domains_idx]
                    ensg_id = fields[gene_idx]
                    transcript_id = fields[transcript_idx]
                    exon_id = fields[exon_idx]
                    
                    if not domains_str or not ensg_id:
                        continue
                    
                    # Process each domain annotation
                    for db_name, hierarchy in self._parse_domain_hierarchy(domains_str):
                        current_dict = self.scores_dict[ensg_id][db_name]
                        
                        # Navigate through domain hierarchy
                        for level in hierarchy:
                            current_dict = current_dict[level]
                            if transcript_id not in current_dict:
                                current_dict[transcript_id] = set()
                            else:
                                current_dict[transcript_id].add(exon_id)
                            # Initialize distribution if not exists
                            if 'distribution' not in current_dict:
                                current_dict['distribution'] = {}
                            # Add score to distribution
                            if prot_var not in current_dict['distribution']:
                                current_dict['distribution'][prot_var] = am_score
                            # elif am_score > current_dict['distribution'][aa_pos]:
                            #     # Keep the highest score for this AA position
                            #     current_dict['distribution'][aa_pos] = am_score
                        
            except (KeyError, ValueError) as e:
                print(f"Error processing record {record.chrom}:{record.pos} - {str(e)}")
                continue

    def finalize_scores(self):
        """Convert score dict to numpy arrays."""
        def convert_nested(d):
            for k, v in d.items():
                if k == 'distribution' and isinstance(v, dict):
                    d[k] = np.array(list(v.values()))
                elif isinstance(v, dict):
                    convert_nested(v)
        
        convert_nested(self.scores_dict)
    
    @staticmethod
    def load_scores_from_pickle(pickle_path: str) -> Dict:
        """Load scores from the pickle file."""
        with open(pickle_path, 'rb') as f:
            data = pickle.load(f)
        return data


def main(vcf_path, output_pickle):
    """
    Notice that the input VCF file must be the AlphaMissense VCF file annotated by VEP with the --domains argument.
    """
    collector = DomainAMScoreCollector(vcf_path)
    collector.collect_scores()
    collector.finalize_scores()
    
    # Output the scores (nested dict) to a pickle file
    if output_pickle:
        with open(output_pickle, 'wb') as f:
            pickle.dump(collector.scores_dict, f)

    return collector.scores_dict
        

if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2])