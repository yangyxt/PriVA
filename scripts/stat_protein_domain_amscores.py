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
    '''
    {
        <ENSG_ID>: {  # Key: Ensembl Gene ID (string, e.g., 'ENSG00000123456')
            <DB_NAME>: {  # Key: Domain Database Name (string, e.g., 'PANTHER', 'Pfam')
            <LEVEL_1_ID>: {  # Key: First level of domain ID (string, e.g., 'PTHR26451', 'PF00001')
                # --- This level contains transcript info and the score distribution ---
                <TRANSCRIPT_ID_1>: set(<EXON_ID_A>, <EXON_ID_B>, ...), # Key: Ensembl Transcript ID (string), Value: Set of Exon IDs (strings like '1/10')
                <TRANSCRIPT_ID_2>: set(<EXON_ID_C>, ...),
                # ... potentially more transcript IDs ...
                'distribution': np.ndarray([<score1>, <score2>, ...]), # Key: literal string 'distribution', Value: NumPy array of AlphaMissense scores (floats) for this domain level

                # --- If the domain has further hierarchy (e.g., PANTHER subfamilies) ---
                <LEVEL_2_ID>: { # Key: Second level of domain ID (string, e.g., 'SF72')
                # --- This level ALSO contains transcript info and its own score distribution ---
                <TRANSCRIPT_ID_1>: set(<EXON_ID_A>, <EXON_ID_D>, ...), # Transcript info specific to this sub-level
                <TRANSCRIPT_ID_3>: set(<EXON_ID_E>, ...),
                # ... potentially more transcript IDs ...
                'distribution': np.ndarray([<score3>, <score4>, ...]), # NumPy array of scores specific to this sub-level (LEVEL_1 + LEVEL_2)

                # --- And potentially more levels ---
                <LEVEL_3_ID>: {
                    # ... structure repeats ...
                    'distribution': np.ndarray([...])
                                }
                            }
                        }
                        # ... potentially more LEVEL_1_IDs ...
                        }
                        # ... potentially more DB_NAMEs ...
                    }
                    # ... potentially more ENSG_IDs ...
    }
    '''
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

                            if 'min_distribution' not in current_dict:
                                current_dict['min_distribution'] = {}
                            if aa_pos not in current_dict['min_distribution']:
                                current_dict['min_distribution'][aa_pos] = am_score
                            elif am_score < current_dict['min_distribution'][aa_pos]:
                                # Keep the lowest score for this AA position
                                current_dict['min_distribution'][aa_pos] = am_score

                            if 'max_distribution' not in current_dict:
                                current_dict['max_distribution'] = {}
                            if aa_pos not in current_dict['max_distribution']:
                                current_dict['max_distribution'][aa_pos] = am_score
                            elif am_score > current_dict['max_distribution'][aa_pos]:
                                # Keep the highest score for this AA position
                                current_dict['max_distribution'][aa_pos] = am_score
            except (KeyError, ValueError) as e:
                print(f"Error processing record {record.chrom}:{record.pos} - {str(e)}")
                continue

    def finalize_scores(self):
        """Convert score dict to numpy arrays."""
        def convert_nested(d):
            for k, v in d.items():
                if k == 'all_distribution' and isinstance(v, dict):
                    d[k] = np.array(list(v.values()))
                elif k == 'min_distribution' and isinstance(v, dict):
                    d[k] = np.array(list(v.values()))
                elif k == 'max_distribution' and isinstance(v, dict):
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