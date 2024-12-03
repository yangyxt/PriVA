#!/usr/bin/env python3

import pysam
from collections import defaultdict
from typing import Dict, List
import sys
import pickle
import logging

logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)
console_handler = logging.StreamHandler()
console_handler.setLevel(logging.INFO)
formatter = logging.Formatter(
    "%(levelname)s:%(asctime)s:%(funcName)s:%(lineno)s:%(message)s"
)
console_handler.setFormatter(formatter)
logger.addHandler(console_handler)

# Define the nested defaultdict factory at the module level
def nested_defaultdict():
    return defaultdict(nested_defaultdict)

class DomainClinVarCollector:
    def __init__(self, vcf_path: str, logger=logger):
        """Initialize collector with VCF file path."""
        self.vcf_path = vcf_path
        self.scores_dict = nested_defaultdict()
        self.logger = logger

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
                hierarchy = parts[1:i + 1]  # Include all levels up to current
                parsed_domains.append((db_name, hierarchy))

        return parsed_domains

    def collect_scores(self):
        """Process VCF file and collect ClinVar annotations per protein domain."""
        vcf = pysam.VariantFile(self.vcf_path)

        # Get CSQ format from header
        csq_format = str(vcf.header).split('Format: ')[-1].strip('>"').split('|')
        domains_idx = csq_format.index('DOMAINS')
        gene_idx = csq_format.index('Gene')
        consq_idx = csq_format.index('Consequence')
        feature_type_idx = csq_format.index('Feature_type')
        transcript_idx = csq_format.index('Feature')

        for record in vcf:
            try:
                if 'CLNREVSTAT' not in record.info or 'CLNSIG' not in record.info:
                    self.logger.warning(
                        f"This record {record.chrom}:{record.pos} is recording an oncogenic variant instead of a Mendelian Pathogenic variant. Continue"
                    )
                    continue

                rev_status = ",".join(record.info['CLNREVSTAT'])
                cln_sig = ",".join(record.info['CLNSIG'])

                csq_value = record.info['CSQ']
                if isinstance(csq_value, tuple):
                    transcript_annotations = csq_value
                else:
                    transcript_annotations = csq_value.split(',')

                for annotation in transcript_annotations:
                    fields = annotation.split('|')

                    # Skip if not a transcript or not an ENST transcript
                    if (
                        fields[feature_type_idx] != 'Transcript'
                        or not fields[transcript_idx].startswith('ENST')
                    ):
                        continue

                    domains_str = fields[domains_idx]
                    # enst_id = fields[transcript_idx]
                    ensg_id = fields[gene_idx]
                    consq = fields[consq_idx]

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
                                current_dict['distribution'] = defaultdict(list)
                            # Collect clinrevstatus, clnsig, and mole_consq
                            current_dict['distribution']["clinrevstatus"].append(rev_status)
                            current_dict['distribution']["clnsig"].append(cln_sig)
                            current_dict['distribution']["mole_consq"].append(consq)

            except (KeyError, ValueError) as e:
                print(f"Error processing record {record.chrom}:{record.pos} - {str(e)}")
                continue

    @staticmethod
    def load_scores_from_pickle(pickle_path: str) -> Dict:
        """Load scores from the pickle file."""
        with open(pickle_path, 'rb') as f:
            data = pickle.load(f)
        return data

def main(vcf_path, output_pickle):
    """
    Notice that the input VCF file must be the ClinVar VCF file annotated by VEP with the --domains argument.
    """
    collector = DomainClinVarCollector(vcf_path)
    collector.collect_scores()

    # Output the scores (nested dict) to a pickle file
    if output_pickle:
        with open(output_pickle, 'wb') as f:
            pickle.dump(collector.scores_dict, f)

    return collector.scores_dict

if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2])