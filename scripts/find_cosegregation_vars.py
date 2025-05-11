import pysam
import pandas as pd
import argparse as ap



def parse_gt(gt):
    """Parse genotype tuple, handling phasing and missing data."""
    if gt is None:
        return None
    alleles = [allele for allele in gt if allele != '.']
    if not alleles:
        return None
    return alleles



def check_variant_presence(alleles, mode, gender, chrom, is_affected):
    """Determine if genotype matches the inheritance mode, considering gender and chromosome."""
    if alleles is None:
        return None
    variant_count = None if alleles is None else sum(1 for allele in alleles if allele == '1')
    is_x = chrom == 'X'
    is_y = chrom == 'Y'

    if is_y:
        if gender == '1':  # Male
            return variant_count >= 1 if is_affected else variant_count == 0
        return False  # Females lack Y chromosome

    if is_x:
        if gender == '1':  # Male, haploid X
            if len(alleles) == 1:
                return alleles[0] == '1' if is_affected else alleles[0] == '0'
            return None
        else:  # Female
            if mode == 'dominant':  # XD
                return variant_count >= 1 if is_affected else variant_count == 0
            elif mode == 'recessive':  # XR
                return all(a == '1' for a in alleles) if is_affected else not all(a == '1' for a in alleles)

    # Autosomal (AD or AR)
    if mode == 'dominant':  # AD
        return variant_count >= 1 if is_affected else variant_count == 0
    elif mode == 'recessive':  # AR
        return all(a == '1' for a in alleles) if is_affected else not all(a == '1' for a in alleles)
    return None



def find_cosegregating_variants(vcf_file, ped_file, mode):
    """Find variants cosegregating under the specified inheritance mode (dominant or recessive).
    
    Args:
        vcf_file (str): Path to multi-family VCF file.
        ped_file (str): Path to pedigree file.
        mode (str): Inheritance mode ('dominant' or 'recessive' or 'both').

    Returns:
        set: Variant tuples (chrom, pos, ref, alt) that cosegregate under the mode.
    """
    if mode not in ['dominant', 'recessive', 'both']:
        raise ValueError("Mode must be 'dominant' or 'recessive' or 'both'")
    modes = ['dominant', 'recessive'] if mode == 'both' else [mode]

    # Load pedigree
    ped_df = pd.read_csv(ped_file, sep='\t', header=None, 
                         names=['FamilyID', 'IndividualID', 'PaternalID', 'MaternalID', 'Sex', 'Phenotype'])
    ped_df.loc[:, 'Sex'] = ped_df['Sex'].astype(int)
    ped_df.loc[:, 'Phenotype'] = ped_df['Phenotype'].astype(int)
    families = ped_df.groupby('FamilyID')

    # Initialize results
    results = {mode: set() for mode in modes}

    # Parse VCF
    with pysam.VariantFile(vcf_file) as vcf:
        for record in vcf:
            chrom = record.chrom
            variant = (chrom, record.pos, record.ref, record.alts[0])
            all_families_pass = True
            for mode in modes:
                for _, family in families:
                    affected = family[family['Phenotype'] == 2][['IndividualID', 'Sex']].values
                    unaffected = family[family['Phenotype'] == 1][['IndividualID', 'Sex']].values

                    # Check affected individuals
                    for ind, gender in affected:
                        gt = record.samples[ind]['GT']
                        alleles = parse_gt(gt)
                        if check_variant_presence(alleles, mode, gender, chrom, True) is False:
                            all_families_pass = False
                            break

                    if not all_families_pass:
                        break

                    # Check unaffected individuals
                    for ind, gender in unaffected:
                        gt = record.samples[ind]['GT']
                        alleles = parse_gt(gt)
                        if check_variant_presence(alleles, mode, gender, chrom, False) is True:
                            all_families_pass = False
                            break

                    if not all_families_pass:
                        break

                if all_families_pass:
                    results[mode].add(variant)
                    
    return results



if __name__ == "__main__":
    parser = ap.ArgumentParser()
    parser.add_argument("--vcf", type=str, required=True)
    parser.add_argument("--ped", type=str, required=True)
    parser.add_argument("--mode", type=str, required=True, choices=['dominant', 'recessive', 'both'])
    args = parser.parse_args()
    results = find_cosegregating_variants(args.vcf, args.ped, args.mode)