import pandas as pd
import argparse as ap
import pysam
import numpy as np
import logging
from enum import Enum
from typing import Optional, Tuple


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
    return len(set(genotype)) == 1 and genotype[0] != None and genotype[0] != 0


def is_variant_carrier(genotype):
    """Check if a sample is a carrier of the variant."""
    return 1 in genotype or 2 in genotype


class VariantSource:
    """
    Class to determine the source of a variant in a proband based on parental genotypes.

    Attributes:
        record (pysam.VariantRecord): The VCF record for the variant.
        proband (str): Sample ID of the proband.
        father (Optional[str]): Sample ID of the father.
        mother (Optional[str]): Sample ID of the mother.
        proband_sex (str): Sex of the proband ('1' for male, '2' for female).
    """

    def __init__(
        self,
        record,
        proband: str,
        father: Optional[str],
        mother: Optional[str],
        proband_sex: str,
    ):
        self.record = record
        self.proband = proband
        self.father = father
        self.mother = mother
        self.proband_sex = proband_sex  # '1' for male, '2' for female
        self.chrom = self.record.chrom
        self.chrom_type = self._get_chrom_type()

    def _get_chrom_type(self) -> str:
        """Determine chromosome type."""
        chrom = self.chrom.replace("chr", "").upper()
        if chrom == "X":
            return "X"
        elif chrom == "Y":
            return "Y"
        else:
            return "AUTOSOMAL"

    def _get_genotype(self, sample: str) -> Optional[Tuple[int, ...]]:
        """Retrieve genotype for a given sample."""
        if sample and sample in self.record.samples:
            gt = self.record.samples[sample].get("GT")
            if gt:
                return gt
        return None

    def _genotype_exists(self, genotype: Optional[Tuple[int, ...]]) -> bool:
        """Check if genotype data exists."""
        return genotype is not None and all(allele is not None for allele in genotype)

    def _has_variant(self, genotype: Optional[Tuple[int, ...]]) -> bool:
        """Check if genotype carries the variant (i.e., has alternate allele)."""
        if not self._genotype_exists(genotype):
            return False
        return 1 in genotype

    def determine_source(self) -> str:
        """
        Determine the variant source.

        Returns:
            str: One of 'PATERNAL', 'MATERNAL', 'DENOVO', 'INHERITED', or 'UNKNOWN'.
        """
        proband_gt = self._get_genotype(self.proband)
        if not self._has_variant(proband_gt):
            return "UNKNOWN"

        father_gt = self._get_genotype(self.father)
        mother_gt = self._get_genotype(self.mother)

        if self.chrom_type == "AUTOSOMAL":
            return self._determine_autosomal_source(proband_gt, father_gt, mother_gt)
        elif self.chrom_type == "X":
            return self._determine_chrX_source(proband_gt, father_gt, mother_gt)
        elif self.chrom_type == "Y":
            return self._determine_chrY_source(proband_gt, father_gt)
        else:
            return "UNKNOWN"

    def _determine_autosomal_source(
        self,
        proband_gt: Tuple[int, ...],
        father_gt: Optional[Tuple[int, ...]],
        mother_gt: Optional[Tuple[int, ...]],
    ) -> str:
        """Determine source for autosomal chromosomes."""
        logger.info(f"For variant {self.record.chrom}:{self.record.pos}:{self.record.id}:{self.record.ref}:{self.record.alts[0]}, proband {self.proband} GT: {proband_gt}, has_variant: {self._has_variant(proband_gt)}, father {self.father} GT: {father_gt}, has_variant: {self._has_variant(father_gt)}, mother {self.mother} GT: {mother_gt}, has_variant: {self._has_variant(mother_gt)}")

        if not self._genotype_exists(father_gt) or not self._genotype_exists(mother_gt):
            # Cannot determine exact source without both parents' genotypes
            if self._has_variant(father_gt) or self._has_variant(mother_gt):
                return "INHERITED"
            else:
                return "UNKNOWN"

        # Both parents' genotypes exist
        proband_has_variant = self._has_variant(proband_gt)
        father_has_variant = self._has_variant(father_gt)
        mother_has_variant = self._has_variant(mother_gt)

        if not father_has_variant and not mother_has_variant:
            return "DENOVO"
        elif father_has_variant and not mother_has_variant:
            return "PATERNAL"
        elif not father_has_variant and mother_has_variant:
            return "MATERNAL"
        elif father_has_variant and mother_has_variant:
            return "INHERITED"
        else:
            return "UNKNOWN"

    def _determine_chrX_source(
        self,
        proband_gt: Tuple[int, ...],
        father_gt: Optional[Tuple[int, ...]],
        mother_gt: Optional[Tuple[int, ...]],
    ) -> str:
        """Determine source for chromosome X."""

        if self.proband_sex == "1":  # Male proband
            # Males have one X chromosome
            if not self._genotype_exists(mother_gt):
                return "UNKNOWN"
            mother_has_variant = self._has_variant(mother_gt)
            if mother_has_variant:
                return "MATERNAL"
            else:
                return "DENOVO"
        elif self.proband_sex == "2":  # Female proband
            return self._determine_autosomal_source(proband_gt, father_gt, mother_gt)
        else:
            return "UNKNOWN"

    def _determine_chrY_source(
        self,
        proband_gt: Tuple[int, ...],
        father_gt: Optional[Tuple[int, ...]],
    ) -> str:
        """Determine source for chromosome Y."""

        if self.proband_sex == "1":  # Male proband
            # Males inherit Y chromosome from father
            if not self._genotype_exists(father_gt):
                return "UNKNOWN"

            father_has_variant = self._has_variant(father_gt)
            if father_has_variant:
                return "PATERNAL"
            else:
                return "DENOVO"
        else:
            # Females do not have Y chromosome
            return "UNKNOWN"


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
    patients = set(fam_ped[fam_ped['Phenotype'] == "2"]['IndividualID'])
    proband = fam_ped[fam_ped['Phenotype'] == "2"]['IndividualID'].tolist()[0]
    proband_sex = fam_ped[fam_ped['IndividualID'] == proband]['Sex'].tolist()[0]
    controls = set(fam_ped[fam_ped['Phenotype'] == "1"]['IndividualID'])
    
    # Define values to be considered as NA
    na_values = ['0', '.', '', 'NA', 'nan', 'NaN', "-", -9]
    ped_df["PaternalID"] = ped_df["PaternalID"].replace(na_values, np.nan)
    ped_df["MaternalID"] = ped_df["MaternalID"].replace(na_values, np.nan)
    # Extract parent IDs
    father = fam_ped["PaternalID"].dropna().drop_duplicates().tolist()[0] if len(fam_ped["PaternalID"].dropna().drop_duplicates().tolist()) > 0 else None
    mother = fam_ped["MaternalID"].dropna().drop_duplicates().tolist()[0] if len(fam_ped["MaternalID"].dropna().drop_duplicates().tolist()) > 0 else None

    # Open input VCF and create output VCF
    with pysam.VariantFile(vcf_path) as vcf_in:
        # Add a variant source INFO field to the header
        vcf_in.header.info.add("VARIANT_SOURCE", "1", "String", "Source of the variant")
        with pysam.VariantFile(output, 'w', header=vcf_in.header) as vcf_out:
            for record in vcf_in:
                # Check if variant should be kept
                keep_variant = True
                if not (1 in record.samples[proband]['GT']):
                    keep_variant = False
                    logger.warning(f"Proband {proband} is not a carrier of the variant {record}")
                    continue

                # Check controls
                for sample in controls:
                    if sample in record.samples:
                        genotype = record.samples[sample]['GT']
                        if is_homozygous_or_hemizygous(genotype):
                            logger.info(f"Control sample {sample} is homozygous or hemizygous for the variant {record}")
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
                        logger.info(f"No patient ({patients}) carries the variant {record}")
                        keep_variant = False
                
                # Write variant to output if it passed filters
                if keep_variant:
                    variant_source = VariantSource(
                        record=record,
                        proband=proband,
                        father=father,
                        mother=mother,
                        proband_sex=str(proband_sex)
                    )
                    source = variant_source.determine_source()
                    # Annotate the source to the INFO field
                    record.info['VARIANT_SOURCE'] = source
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
    parser.add_argument("-o", "--output", type=str, help="The path to the output VCF", required=True)
    
    args=parser.parse_args()
    pedigree_filter(vcf_path=args.vcf,
                    ped_path=args.ped_table,
                    target_fam=args.target_fam,
                    output=args.output)
    
    
