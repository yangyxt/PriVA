#!/usr/bin/env python3
"""
Clean, functional script for determining cis/trans relationships between variants.
Uses gene-specific parallelization with simple, picklable data structures.
"""

import logging
from dataclasses import dataclass
from typing import Optional, List, Dict, Any, Tuple, Set
from collections import defaultdict
import pandas as pd
import numpy as np
import multiprocessing as mp

logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)
console_handler = logging.StreamHandler()
console_handler.setLevel(logging.INFO)
formatter = logging.Formatter("%(levelname)s:%(asctime)s:%(funcName)s:%(lineno)s:%(message)s")
console_handler.setFormatter(formatter)
logger.addHandler(console_handler)


@dataclass
class Genotype:
    """Represents a parsed genotype with phasing information."""
    allele1: Optional[str]
    allele2: Optional[str]
    is_phased: bool
    is_hemizygous: bool = False
    chromosome: str = None
    
    @classmethod
    def parse(cls, gt_str: Optional[str], patient_gender: str, chromosome: str) -> 'Genotype':
        """Parse genotype string into Genotype object."""
        if not gt_str or gt_str in ['.', './.', '.|.']:
            return cls(None, None, False, False, chromosome)
        
        delimiter = gt_str[1]
        
        # Handle haploid genotypes (hemizygous)
        if (patient_gender == 'male' and ('X' in chromosome)) or ("Y" in chromosome) or ("M" in chromosome):
            alleles = gt_str.split(delimiter)

            if len(alleles) == 1 and alleles[0] != '.':
                return cls(alleles[0], alleles[0], True, True, chromosome)
            
            allele1, allele2 = alleles

            # Handle missing alleles
            if allele1 == '.':
                return cls(allele2, allele2, True, True, chromosome)
            elif allele2 == '.':
                return cls(allele1, allele1, True, True, chromosome)
            elif allele1 == allele2:
                return cls(allele1, allele2, True, True, chromosome)
            else:
                if "1" in alleles:
                    return cls("1", "1", True, True, chromosome)
                raise ValueError(f"Invalid genotype: {gt_str} for variant in patient (gender {patient_gender}) at chromosome {chromosome}")
            
        # Handle diploid genotypes
        if delimiter == '|':
            is_phased = True
        else:
            is_phased = False

        alleles = gt_str.split(delimiter)
        if len(alleles) != 2:
            raise ValueError(f"Invalid genotype: {gt_str} for variant in patient (gender {patient_gender}) at chromosome {chromosome}")
            
        allele1 = alleles[0] if alleles[0] != '.' else None
        allele2 = alleles[1] if alleles[1] != '.' else None
        
        return cls(allele1, allele2, is_phased, False, chromosome)
    
    
    @property
    def has_alt(self) -> bool:
        """Check if genotype has alternative allele."""
        return self.allele1 == '1' or self.allele2 == '1'
    
    @property
    def full_genotype(self) -> bool:
        """Check if genotype is valid for analysis."""
        if self.is_hemizygous:
            return (self.allele1 is not None) or (self.allele2 is not None)
        else:
            return (self.allele1 is not None) and (self.allele2 is not None)


@dataclass
class Variant:
    """Represents a genomic variant with genotype information."""
    variant_id: str
    chromosome: str
    position: int
    ref: str
    alt: str
    gene: str
    is_pathogenic: bool
    patient_id: str
    patient_gender: str
    patient_gt: Genotype
    father_gt: Optional[Genotype] = None
    mother_gt: Optional[Genotype] = None
    
    def __post_init__(self):
        """Normalize chromosome and gender."""
        self.chromosome = str(self.chromosome)
        self.patient_gender = "male" if self.patient_gender == "1" else "female"
    
    def __hash__(self) -> int:
        """Make Variant hashable by variant_id for use in sets/dicts."""
        return hash(self.variant_id)
    
    
    def determine_phase_with(self, other: 'Variant') -> str:
        """
        Determine phase relationship with another variant.
        
        Returns: 'cis', 'trans', or 'unknown', or 'both' for homozygous variants
         """
        # Basic validation
        if (self.patient_id != other.patient_id or 
            self.chromosome != other.chromosome or
            not self.patient_gt.has_alt or 
            not other.patient_gt.has_alt):
            return 'unknown'
        
        # Male X variants are always in cis (hemizygous)
        if self.patient_gt.is_hemizygous and other.patient_gt.is_hemizygous and \
           self.patient_gt.has_alt and other.patient_gt.has_alt:
            return 'cis'
        
        # If both are phased, determine directly
        if self.patient_gt.is_phased and other.patient_gt.is_phased:
            return self._phase_from_haplotypes(other)
        
        # Use parental information for unphased variants
        return self._phase_with_parents(other)
    
    
    def _phase_from_haplotypes(self, other: 'Variant') -> str:
        """Determine phase from phased genotypes."""
        # Deal with homozygous case
        self_hom = (self.patient_gt.allele1 == '1' and self.patient_gt.allele2 == '1')
        other_hom = (other.patient_gt.allele1 == '1' and other.patient_gt.allele2 == '1')
        
        if self_hom or other_hom:
            return 'both'
        
        # Now both are heterozygous variants
        if (self.patient_gt.allele1 == other.patient_gt.allele1) or (self.patient_gt.allele2 == other.patient_gt.allele2):
            return "cis"
        elif (self.patient_gt.allele1 == other.patient_gt.allele2) or (self.patient_gt.allele2 == other.patient_gt.allele1):
            return "trans"
        else:
            raise ValueError(f"Invalid genotype: {self.patient_gt} for variant in patient (gender {self.patient_gender}) at chromosome {self.chromosome}, please try to convert the genotype to bi-alleleic first")
        
    
    def _phase_with_parents(self, other: 'Variant') -> str:
        """Use parental genotypes to determine phase."""
        # Check if we have complete parental information for both variants
        if (self.father_gt is None or self.mother_gt is None or 
            other.father_gt is None or other.mother_gt is None):
            return 'unknown'
        
        # Check if parents have complete genotypes
        if (not self.father_gt.full_genotype or not self.mother_gt.full_genotype or
            not other.father_gt.full_genotype or not other.mother_gt.full_genotype):
            return 'unknown'
        
        # Determine which parents could have transmitted the alt allele for each variant
        self_from_father = self.father_gt.has_alt
        self_from_mother = self.mother_gt.has_alt
        other_from_father = other.father_gt.has_alt
        other_from_mother = other.mother_gt.has_alt
        
        # If neither parent has the alt allele, something is wrong
        if not self_from_father and not self_from_mother:
            return 'unknown'
        if not other_from_father and not other_from_mother:
            return 'unknown'
        
        # Case 1: self variant must come from father only
        if self_from_father and not self_from_mother:
            if other_from_father and not other_from_mother:
                # Both variants must come from father → cis
                return 'cis'
            elif not other_from_father and other_from_mother:
                # self from father, other from mother → trans
                return 'trans'
            else:
                # other variant could come from either parent → unknown
                return 'unknown'
        
        # Case 2: self variant must come from mother only
        elif not self_from_father and self_from_mother:
            if not other_from_father and other_from_mother:
                # Both variants must come from mother → cis
                return 'cis'
            elif other_from_father and not other_from_mother:
                # self from mother, other from father → trans
                return 'trans'
            else:
                # other variant could come from either parent → unknown
                return 'unknown'
        
        # Case 3: self variant could come from either parent
        else:  # self_from_father and self_from_mother
            if other_from_father and not other_from_mother:
                # other must come from father only
                # If self is heterozygous and other must be from father,
                # we can't determine phase without additional information
                return 'unknown'
            elif not other_from_father and other_from_mother:
                # other must come from mother only
                # If self is heterozygous and other must be from mother,
                # we can't determine phase without additional information
                return 'unknown'
            else:
                # Both variants could come from either parent → unknown
                return 'unknown'


def create_variants_from_gene_df(gene_df: pd.DataFrame, gene_symbol: str, patient_info: Dict[str, Dict[str, Any]]) -> List[Variant]:
    """
    Create Variant objects from a gene-specific DataFrame.
    
    Args:
        gene_df: DataFrame containing all rows for a specific gene
        gene_symbol: Gene symbol being processed
        patient_info: Dictionary mapping patient_id -> patient information
        
    Returns:
        List of Variant objects for this gene
    """
    patient_ids = set(patient_info.keys())
    patient_variants = {pid: {"pathogenic": [], "non_pathogenic": []} for pid in patient_ids}
    
    # Group by variant_id to get unique genomic variants
    for variant_id, variant_group in gene_df.groupby('variant_id'):
        if pd.isna(variant_id):
            continue
        
        # Check if any transcript is pathogenic
        is_pathogenic = variant_group['is_pathogenic'].any()
        
        # Get variant genomic information from first row
        first_row = variant_group.iloc[0]
        chromosome = str(first_row['chrom'])
        position = int(first_row['pos'])
        ref = str(first_row['ref'])
        alt = str(first_row['alt'])
        
        # Create variant for each patient that has this variant
        for patient_id, patient_data in patient_info.items():
            if patient_id not in gene_df.columns:
                continue
                
            # Get patient's genotype for this variant
            patient_gt_str = first_row.get(patient_id)
            if pd.isna(patient_gt_str):
                continue
            
            # Parse genotypes
            patient_gt = Genotype.parse(patient_gt_str, patient_data['gender'], chromosome)
            
            # Skip if no alternative allele or invalid genotype
            if not patient_gt.has_alt:
                continue
            
            # Parse parental genotypes if available
            father_gt = None
            mother_gt = None
            
            if patient_data.get('father_id') and patient_data['father_id'] in gene_df.columns:
                father_gt_str = first_row.get(patient_data['father_id'])
                if not pd.isna(father_gt_str):
                    father_gt = Genotype.parse(father_gt_str, "male", chromosome)
            
            if patient_data.get('mother_id') and patient_data['mother_id'] in gene_df.columns:
                mother_gt_str = first_row.get(patient_data['mother_id'])
                if not pd.isna(mother_gt_str):
                    mother_gt = Genotype.parse(mother_gt_str, "female", chromosome)  # Mother is always female
            
            # Create variant object
            variant = Variant(
                variant_id=str(variant_id),
                chromosome=chromosome,
                position=position,
                ref=ref,
                alt=alt,
                gene=gene_symbol,
                is_pathogenic=is_pathogenic,
                patient_id=patient_id,
                patient_gender=patient_data['gender'],
                patient_gt=patient_gt,
                father_gt=father_gt,
                mother_gt=mother_gt
            )

            if is_pathogenic:
                patient_variants[patient_id]["pathogenic"].append(variant)
            else:
                patient_variants[patient_id]["non_pathogenic"].append(variant)
    
    return patient_variants


def analyze_gene_variants(gene_df: pd.DataFrame, gene_symbol: str, patient_info: Dict[str, Dict[str, Any]]) -> Tuple[Dict[str, Set[str]], Dict[str, Set[str]]]:
    """
    Analyze variants within a gene to find cis/trans relationships.
    Uses optimized algorithm: each variant pairs with all pathogenic variants.
    
    Args:
        gene_df: DataFrame for a specific gene
        gene_symbol: Gene symbol being analyzed
        patient_info: Patient information dictionary
        
    Returns:
        Tuple of (cis_variants_per_patient, trans_variants_per_patient)
        where each dict maps patient_id -> Set[variant_ids]
    """
    # Create all variant objects for this gene
    variants_by_patient = create_variants_from_gene_df(gene_df, gene_symbol, patient_info)

    # Track results per patient
    variants_phase_patho = defaultdict(dict)
    
    # Track processed pairs to avoid redundant comparisons
    processed_pairs = defaultdict(dict)
    
    # For each patient, analyze their variants
    for patient_id, patient_variants in variants_by_patient.items():
        if len(patient_variants["pathogenic"]) == 0:
            continue
        
        # For each variant, compare with all pathogenic variants
        for pvar in patient_variants["pathogenic"]:
            for upvar in patient_variants["non_pathogenic"]:
                # Determine phase relationship
                phase = None
                if pvar in processed_pairs:
                    if upvar in processed_pairs[pvar]:
                        phase = processed_pairs[pvar][upvar]
                
                if phase is None:
                    phase = upvar.determine_phase_with(pvar)
                    processed_pairs[pvar.variant_id][upvar.variant_id] = phase
                    processed_pairs[upvar.variant_id][pvar.variant_id] = phase
                
                variants_phase_patho[patient_id][upvar.variant_id] = phase
    
    return variants_phase_patho


def process_gene_wrapper(args: Tuple[str, pd.DataFrame, Dict[str, Dict[str, Any]]]) -> Tuple[str, Dict[str, Set[str]], Dict[str, Set[str]]]:
    """
    Wrapper function for multiprocessing gene analysis.
    
    Args:
        args: Tuple of (gene_symbol, gene_df, patient_info)
        
    Returns:
        Tuple of (gene_symbol, cis_variants_per_patient, trans_variants_per_patient)
    """
    gene_symbol, gene_df, patient_info = args
    variants_phase_patho = analyze_gene_variants(gene_df, gene_symbol, patient_info)
    return variants_phase_patho


def build_patient_info(pedigree_df: pd.DataFrame, variants_df: pd.DataFrame) -> Dict[str, Dict[str, Any]]:
    """Build patient information mapping."""
    sample_ids = set(variants_df.columns.tolist())
    pedigree_df = pedigree_df.loc[pedigree_df['IndividualID'].isin(sample_ids), :]
    patients = pedigree_df.loc[pedigree_df["Phenotype"].isin(["2", 2]), "IndividualID"].unique().tolist()
    
    patient_info = {}
    for patient_id in patients:
        row = pedigree_df.loc[pedigree_df["IndividualID"] == patient_id, :].squeeze()

        # Normalize gender
        sex_val = row.get('Sex')
        gender = "male" if sex_val in [1, "1"] else "female" if sex_val in [2, "2"] else None
        if gender is None:
            logger.error(f"Unknown gender for patient from pedigree table {patient_id}: {sex_val}")
            # Check whether variants_df contain chrY
            if "chrY" in variants_df.loc[variants_df[patient_id].str.contains("1"), 'chrom'].unique().tolist():
                gender = "male"
            else:
                gender = "female"
            logger.warning(f"Set gender for patient {patient_id} to {gender} based on variants_df")
        
        # Get parent IDs
        father_id = row.get('PaternalID')
        mother_id = row.get('MaternalID')
        
        father_id = str(father_id) if pd.notna(father_id) and father_id not in [0, '0'] else None
        mother_id = str(mother_id) if pd.notna(mother_id) and mother_id not in [0, '0'] else None
        
        patient_info[patient_id] = {
            'gender': gender,
            'father_id': father_id,
            'mother_id': mother_id,
            'phenotype': row.get('Phenotype')
        }
    
    return patient_info


def determine_cis_trans_relationships(
    variants_df: pd.DataFrame,
    pathogenic_array: List[bool],
    pedigree_df: pd.DataFrame,
    num_processes: Optional[int] = None
) -> Tuple[np.ndarray, np.ndarray]:
    """
    Main function to determine cis/trans relationships using gene-specific parallelization.
    
    Args:
        variants_df: DataFrame with variant annotations  
        pathogenic_array: Boolean array for pathogenic status
        pedigree_df: Pedigree information
        patient_phenotype_value: Phenotype value for patients
        num_processes: Number of processes for multiprocessing
        
    Returns:
        Tuple of (cis_count_array, trans_count_array)
    """
    # Prepare data
    variants_df['is_pathogenic'] = pathogenic_array
    
    # Build patient info
    patient_info = build_patient_info(pedigree_df, variants_df)
    
    logger.info(f"Processing {len(patient_info)} patients across genes")
    
    uniq_genes = variants_df['Gene'].unique().tolist()
    logger.info(f"Processing {len(uniq_genes)} genes")

    arg_generator = ((gene_symbol, variants_df.loc[variants_df['Gene'] == gene_symbol, :], patient_info) for gene_symbol in uniq_genes)
    
    # Process genes (with optional multiprocessing)
    pool = None
    if num_processes and num_processes > 1 and len(uniq_genes) > 1:
        pool = mp.Pool(num_processes)
        results = pool.map(process_gene_wrapper, arg_generator)
    else:
        results = list(map(process_gene_wrapper, arg_generator))

    total_phase_map = defaultdict(dict)
    for variants_phase_patho in results:
        for patient_id in patient_info.keys():
            total_phase_map[patient_id] = { **total_phase_map[patient_id], **variants_phase_patho[patient_id] }
    
    if pool is not None:
        pool.close()
        pool.join()

    # Aggregate results from all genes
    in_cis_pathogenic = np.zeros(len(variants_df), dtype=int)
    in_trans_pathogenic = np.zeros(len(variants_df), dtype=int)
    for patient_id in patient_info.keys():
        variants_df[patient_id + "_phase_patho"] = variants_df["variant_id"].map(total_phase_map[patient_id])
        in_cis_pathogenic = np.where(variants_df[patient_id + "_phase_patho"].isin(["cis", "both"]).to_numpy(), in_cis_pathogenic + 1, in_cis_pathogenic)
        in_trans_pathogenic = np.where(variants_df[patient_id + "_phase_patho"].isin(["trans", "both"]).to_numpy(), in_trans_pathogenic + 1, in_trans_pathogenic)

    return in_cis_pathogenic, in_trans_pathogenic, variants_df

