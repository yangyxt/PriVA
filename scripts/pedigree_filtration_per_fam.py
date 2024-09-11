import pandas as pd
import argparse as ap
import pysam
import os
import sys
import numpy as np
import re
import logging


logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)
console_handler=logging.StreamHandler()
console_handler.setLevel(logging.INFO)
formatter = logging.Formatter("%(levelname)s:%(asctime)s:%(funcName)s:%(lineno)s:%(message)s")
console_handler.setFormatter(formatter)
logger.addHandler(console_handler)


# Establish a modularized script to filter the variants on a VCF file (only contain variants from samples within the same family)
# Use pedigree table (can be read by pandas) (the columns are #FamilyID IndividualID PaternalID MaternalID Sex (1 for male and 2 for female) Phenotype) to decide which samples are healthy (Phenotype value is 1) and which samples are ill (Phenotype value is 2)
# Use pysam to parse the vcf file and perform pedigree based filtration on each variant record
# Each filtration should be executed by a centralized function (returning boolean values) designed for this single purpose. 
# By principle, we only filter out the variant records that:
# 1. homozygous or hemizygous in any control samples
# 2. not carried by any patients in this family

def is_homozygous_or_hemizygous(genotype):
    """Check if a genotype is homozygous or hemizygous."""
    alleles = genotype.split(genotype[1])
    return len(alleles) == 1 and alleles[0] != '.' and alleles[0] != "0"

def is_variant_carrier(genotype):
    """Check if a sample is a carrier of the variant."""
    return '1' in genotype or '2' in genotype

def pedigree_filter(vcf_path, ped_path, target_fam, output):
    """
    Filter variants based on pedigree information.
    
    Args:
    vcf_path (str): Path to the input VCF file.
    ped_path (str): Path to the pedigree file.
    target_fam (str): Target family ID.
    output (str): Path to the output VCF file.
    """
    # Read pedigree file
    ped_df = pd.read_csv(ped_path, sep='\t', header=None, 
                         names=['#FamilyID', 'IndividualID', 'PaternalID', 'MaternalID', 'Sex', 'Phenotype'])
    
    # Filter for target family
    fam_ped = ped_df[ped_df['#FamilyID'] == target_fam]
    
    if fam_ped.empty:
        logger.error(f"Family {target_fam} not found in pedigree file.")
        return
    
    # Identify patients and controls
    patients = set(fam_ped[fam_ped['Phenotype'] == 2]['IndividualID'])
    controls = set(fam_ped[fam_ped['Phenotype'] == 1]['IndividualID'])
    
    # Open input VCF and create output VCF
    with pysam.VariantFile(vcf_path) as vcf_in, pysam.VariantFile(output, 'w', header=vcf_in.header) as vcf_out:
        for record in vcf_in:
            # Check if variant should be kept
            keep_variant = True
            
            # Check controls
            for sample in controls:
                if sample in record.samples:
                    genotype = record.samples[sample]['GT']
                    if is_homozygous_or_hemizygous(genotype):
                        keep_variant = False
                        break
            
            # Check patients
            if keep_variant:
                patient_carrier = False
                for sample in patients:
                    if sample in record.samples:
                        genotype = record.samples[sample]['GT']
                        if is_variant_carrier(genotype):
                            patient_carrier = True
                            break
                
                if not patient_carrier:
                    keep_variant = False
            
            # Write variant to output if it passed filters
            if keep_variant:
                vcf_out.write(record)
    
    logger.info(f"Filtered VCF written to {output}")



if __name__ == '__main__':
    # args.var_table should be the path to the table of variants.
    # args.ped_table should be the path to the table of pedigree info.
    # args.target_fam should be the name of family ID.
    # args.output should be the path to the output table.
    logger = logging.getLogger(__name__)
    logger.setLevel(logging.INFO)
    
    parser = ap.ArgumentParser()
    parser.add_argument("-v", "--vcf", type=str, help="The path to the VCF file", required=True)
    parser.add_argument("-p", "--ped_table", type=str, help="The path to the table of pedigree info", required=True)
    parser.add_argument("-f", "--target_fam", type=str, help="The Family ID of the target fam in this batch", required=True)
    parser.add_argument("-o", "--output", type=str, help="The path to the output table", required=True)
    
    args=parser.parse_args()
    pedigree_filter(var_path=args.var_table,ped_path=args.ped_table,target_fam=args.target_fam,output=args.output)
    
    
