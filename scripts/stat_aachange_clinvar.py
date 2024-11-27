#!/usr/bin/env python3

import pysam
from collections import defaultdict
from typing import Dict
import sys
import pickle
import logging

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
        self.logger = logger

    def collect_clinvar_aa(self):
        """Process VCF file and collect ClinVar amino acid changes."""
        vcf = pysam.VariantFile(self.vcf_path)

        # Get CSQ format from header
        csq_format = str(vcf.header).split('Format: ')[-1].strip('>"').split('|')
        aa_pos_idx = csq_format.index('Protein_position')
        hgvsp_idx = csq_format.index('HGVSp')
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

                    # Skip if not a transcript
                    if fields[feature_type_idx] != 'Transcript':
                        continue

                    enst_id = fields[transcript_idx]
                    consq = fields[consq_idx]
                    aa_pos = fields[aa_pos_idx]
                    hgvsp = fields[hgvsp_idx]

                    if not aa_pos or not hgvsp or not enst_id:
                        continue

                    current_dict = self.clinvaa_dict[enst_id][aa_pos][hgvsp]

                    if "CLNSIG" not in current_dict:
                        current_dict["CLNSIG"] = [cln_sig]
                    else:
                        current_dict["CLNSIG"].append(cln_sig)

                    if "CLNREVSTAT" not in current_dict:
                        current_dict["CLNREVSTAT"] = [rev_status]
                    else:
                        current_dict["CLNREVSTAT"].append(rev_status)

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
    Notice that the input VCF file must be the ClinVar VCF file annotated by VEP.
    """
    collector = AAClinvarCollector(vcf_path)
    collector.collect_clinvar_aa()

    # Output the scores (nested dict) to a pickle file
    if output_pickle:
        with open(output_pickle, 'wb') as f:
            pickle.dump(collector.clinvaa_dict, f)

    return collector.clinvaa_dict

if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2])