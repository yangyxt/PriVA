#!/usr/bin/env python3

import pysam
import numpy as np
from collections import defaultdict
from typing import Dict, List
import sys
import json


class DomainScoreCollector:
    def __init__(self, vcf_path: str):
        """Initialize collector with VCF file path."""
        self.vcf_path = vcf_path
        self.scores_dict = self._create_nested_dict()
        
    @staticmethod
    def _create_nested_dict():
        """Create a nested defaultdict that allows arbitrary depth."""
        return defaultdict(lambda: defaultdict(DomainScoreCollector._create_nested_dict))
    
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
        feature_type_idx = csq_format.index('Feature_type')
        transcript_idx = csq_format.index('Feature')
        
        for record in vcf:
            try:
                am_score = float(record.info['AM_PATHOGENICITY'])
                prot_var = record.info['PVAR']
                aa_pos = prot_var[:-1]
                
                # Handle CSQ as tuple - take first element if it's a tuple
                csq_value = record.info['CSQ']
                if isinstance(csq_value, tuple):
                    transcript_annotations = csq_value
                else:
                    transcript_annotations = csq_value.split(',')
                
                for annotation in transcript_annotations:
                    fields = annotation.split('|')
                    
                    # Skip if not a transcript or not an ENST transcript
                    if (fields[feature_type_idx] != 'Transcript' or 
                        not fields[transcript_idx].startswith('ENST')):
                        continue
                    
                    domains_str = fields[domains_idx]
                    ensg_id = fields[gene_idx]
                    
                    if not domains_str or not ensg_id:
                        continue
                    
                    # Process each domain annotation
                    for db_name, hierarchy in self._parse_domain_hierarchy(domains_str):
                        current_dict = self.scores_dict[ensg_id][db_name]
                        
                        # Navigate through domain hierarchy
                        for level in hierarchy:
                            current_dict = current_dict[level]
                            # Initialize distribution if not exists
                            if 'distribution' not in current_dict:
                                current_dict['distribution'] = {}
                            # Add score to distribution
                            if aa_pos not in current_dict['distribution']:
                                current_dict['distribution'][aa_pos] = am_score
                            elif am_score > current_dict['distribution'][aa_pos]:
                                # We only keep the strongest structural change at this AA position
                                # The resulting distribution will follow Poisson distribution
                                current_dict['distribution'][aa_pos] = am_score
                        
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
    
    def convert_for_json(self, d) -> Dict:
        """Return the collected scores dictionary with numpy arrays converted to lists."""
        output = {}
        for k, v in d.items():
            if isinstance(v, dict):
                output[k] = self.convert_for_json(v)
            elif isinstance(v, np.ndarray):
                output[k] = v.tolist()  # Convert numpy array to list
            else:
                output[k] = v
        return output
    

    def load_scores_from_json(json_path: str) -> Dict:
        """Load scores from JSON and convert lists back to numpy arrays."""
        def convert_from_json(d):
            output = {}
            for k, v in d.items():
                if isinstance(v, dict):
                    output[k] = convert_from_json(v)
                elif k == 'distribution' and isinstance(v, list):
                    output[k] = np.array(v)  # Convert list back to numpy array
                else:
                    output[k] = v
            return output
        
        with open(json_path, 'r') as f:
            data = json.load(f)
        return convert_from_json(data)


def main(vcf_path, output_json):
    '''
    Notice that the input VCF file must be the AlphaMissense VCF file annotated by VEP with the --domains argument.
    '''
    collector = DomainScoreCollector(vcf_path)
    collector.collect_scores()
    collector.finalize_scores()
    
    # Get scores with numpy arrays converted to lists
    scores_dict = collector.convert_for_json(collector.scores_dict)

    # Output the scores (nested dict) to a JSON file
    if output_json:
        with open(output_json, 'w') as f:
            json.dump(scores_dict, f, indent=4)

    return collector.scores_dict
        
    

if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2])