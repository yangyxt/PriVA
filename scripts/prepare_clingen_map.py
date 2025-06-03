import pandas as pd
import pysam
import numpy as np
import gzip
import pickle
import argparse
import sys

def create_clinvar_maps(vcf_path):
    """
    Create dictionaries mapping ClinVar IDs to variant IDs
    
    Args:
        vcf_path: Path to the ClinVar VCF file
        
    Returns:
        Dictionary mapping ClinVar ID to variant_id
    """
    # Initialize dictionary - reversed mapping
    clinvar_id_map = {}   # ClinVar ID -> variant_id
    
    # Open VCF file
    vcf = pysam.VariantFile(vcf_path)
    
    # Process each record
    for record in vcf:
        # Create variant ID
        variant_id = f"{record.chrom}:{record.pos}:{record.ref}-{record.alts[0]}"
        
        # Handle record.id - it can be None or a string
        if record.id is not None:
            # record.id is a string, could be multiple IDs separated by semicolons
            clinvar_ids = str(record.id).split(';')
            for clinvar_id in clinvar_ids:
                clinvar_id = clinvar_id.strip()
                if clinvar_id and clinvar_id != '.':  # Skip empty or missing IDs
                    clinvar_id_map[clinvar_id] = variant_id
    
    return clinvar_id_map


def prepare_clingen_map(clinvar_vcf, clingen_tab, output_pkl):
	clingen_df = pd.read_table(clingen_tab, low_memory=False)
	clingen_df.dropna(subset=['ClinVar Variation Id'], inplace=True)
	clingen_df = clingen_df[['ClinVar Variation Id', 'Applied Evidence Codes (Met)']]
	clingen_df.loc[:, "ClinVar Variation Id"] = clingen_df["ClinVar Variation Id"].astype(int).astype(str)
	print(clingen_df.head().to_string(index=False), file=sys.stderr)
	clinvar_id_map = create_clinvar_maps(clinvar_vcf)
	clingen_df["variant_id"] = clingen_df["ClinVar Variation Id"].map(clinvar_id_map)
	clingen_df.dropna(subset=['variant_id'], inplace=True)
	print(clingen_df.head().to_string(index=False), file=sys.stderr)
	varid_acmg_map = clingen_df.set_index("variant_id")['Applied Evidence Codes (Met)'].to_dict()
	with gzip.open(output_pkl, 'wb') as f:
		pickle.dump(varid_acmg_map, f)


if __name__ == "__main__":
	parser = argparse.ArgumentParser()
	parser.add_argument("-v", "--clinvar_vcf", type=str, required=True)
	parser.add_argument("-t", "--clingen_tab", type=str, required=True)
	parser.add_argument("-o", "--output_pkl", type=str, required=True)
	args = parser.parse_args()
	prepare_clingen_map(args.clinvar_vcf, args.clingen_tab, args.output_pkl)

