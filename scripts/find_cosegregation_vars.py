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



def check_variant_presence(alleles, gender, chrom):
    """Determine if genotype matches the inheritance mode, considering gender and chromosome."""
    if alleles is None:
        return None
    variant_count = None if alleles is None else sum(1 for allele in alleles if allele == '1')
    is_x = 'X' in chrom
    is_y = 'Y' in chrom

    if is_y:
        if gender == '1':  # Male
            return variant_count >= 1
        return False  # Females lack Y chromosome

    if is_x:
        if gender == '1':  # Male, haploid X
            return "Hemizygous" if variant_count == 1 else False
        else:  # Female
            if variant_count > 1:
                return "Homozygous"
            elif variant_count == 1:
                return "Heterozygous"
            else:
                return False

    # Autosomal (AD or AR)
    if variant_count > 1:
        return "Homozygous"
    elif variant_count == 1:
        return "Heterozygous"
    else:
        return False



def find_cosegregating_variants(vcf_file, ped_file):
    """Find variants cosegregating under the specified inheritance mode (dominant or recessive).
    
    Args:
        vcf_file (str): Path to multi-family VCF file.
        ped_file (str): Path to pedigree file.
    Returns:
        set: Variant tuples (chrom, pos, ref, alt) that cosegregate under the mode.
    """

    # Load pedigree
    ped_df = pd.read_csv(ped_file, sep='\t', header=None, 
                         names=['FamilyID', 'IndividualID', 'PaternalID', 'MaternalID', 'Sex', 'Phenotype'])
    ped_df.loc[:, 'Sex'] = ped_df['Sex'].astype(int)
    ped_df.loc[:, 'Phenotype'] = ped_df['Phenotype'].astype(int)
    families = ped_df.groupby('FamilyID')

    # Per-var stats
    per_var_stats = {}

    # Parse VCF
    with pysam.VariantFile(vcf_file) as vcf:
        for record in vcf:
            chrom = record.chrom
            variant = (chrom, record.pos, record.ref, record.alts[0])
            all_families_pass = True
            affected_segregation_count = 0
            unaffected_segregation_count = 0
            affected_homo_count = 0
            unaffected_nonhomo_count = 0
            # Initialize results
            results = { "Affected_segregated_inds": 0,
                        "Unaffected_segregated_inds": 0,
                        "Affected_homozygous_inds": 0,
                        "Unaffected_nonhomo_inds": 0 }
    
            for _, family in families:
                affected = family[family['Phenotype'] == 2][['IndividualID', 'Sex']].values
                unaffected = family[family['Phenotype'] == 1][['IndividualID', 'Sex']].values

                # Identify parent sample ID (if available)
                father = [x for x in family['PaternalID'].dropna().drop_duplicates().values if x not in [0, "0", ".", "-"]]
                father = father[0] if father else None
                mother = [x for x in family['MaternalID'].dropna().drop_duplicates().values if x not in [0, "0", ".", "-"]]
                mother = mother[0] if mother else None
                parents = []
                parents = parents + [father] if father else parents
                parents = parents + [mother] if mother else parents

                # Check affected individuals
                affected_ind_gts = []
                male_affected_gts = []
                for ind, gender in affected:
                    gt = record.samples[ind]['GT']
                    alleles = parse_gt(gt)
                    presence = check_variant_presence(alleles, gender, chrom)
                    if not presence:
                        all_families_pass = False
                        break
                    else:
                        affected_ind_gts.append(presence)
                        if gender == 1:
                            male_affected_gts.append(presence)

                if not all_families_pass:
                    # Early exit if affected individuals do not carry the variant
                    break

                # Check unaffected individuals
                non_parental_controls_gts = []
                parental_controls_gts = []
                male_nonparental_control_gts = []

                for ind, gender in unaffected:
                    gt = record.samples[ind]['GT']
                    alleles = parse_gt(gt)
                    presence = check_variant_presence(alleles, gender, chrom)
                    if presence in [ "Homozygous", "Hemizygous"]:
                        all_families_pass = False
                        break
                    else:
                        if ind in parents:
                            parental_controls_gts.append(presence)
                        else:
                            non_parental_controls_gts.append(presence)
                            if gender == 1:
                                male_nonparental_control_gts.append(presence)

                if not all_families_pass:
                    # Early exit if any unaffected individual carries the variant in disallowed state
                    break

                control_gts = parental_controls_gts + non_parental_controls_gts

                if not any(control_gts):
                    affected_segregation_count += len(affected_ind_gts)
                    unaffected_segregation_count += len(non_parental_controls_gts)
                    affected_male_segregation_count += len(male_affected_gts)
                    unaffected_male_sgregation_count += len(male_nonparental_control_gts)
                elif "Heterozygous" in affected_ind_gts:
                    all_families_pass = False
                else:
                    affected_homo_count += len(affected_ind_gts)
                    unaffected_nonhomo_count += len(non_parental_controls_gts)
                if not all_families_pass:
                    break

            if all_families_pass:
                results["Affected_segregated_inds"] += affected_segregation_count
                results["Affected_homozygous_inds"] += affected_homo_count
                results["Unaffected_segregated_inds"] += unaffected_segregation_count
                results["Unaffected_nonhomo_inds"] += unaffected_nonhomo_count
                results["Male_affected_segregated_inds"] += affected_male_segregation_count
                results["Male_unaffected_segregated_inds"] += unaffected_male_sgregation_count
                per_var_stats[variant] = results
               
    return per_var_stats



if __name__ == "__main__":
    parser = ap.ArgumentParser()
    parser.add_argument("--vcf", type=str, required=True)
    parser.add_argument("--ped", type=str, required=True)
    parser.add_argument("--mode", type=str, required=True, choices=['dominant', 'recessive', 'both'])
    args = parser.parse_args()
    results = find_cosegregating_variants(args.vcf, args.ped, args.mode)