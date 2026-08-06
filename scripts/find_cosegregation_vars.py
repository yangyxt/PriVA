import pysam
import pandas as pd
import argparse as ap



def parse_gt(gt):
    """Parse genotype tuple, handling phasing and missing data."""
    if gt is None:
        return None
    alleles = [allele for allele in gt if allele not in {None, '.'}]
    if not alleles:
        return None
    return alleles


def _normalized_sex(gender):
    value = str(gender).strip().lower()
    if value in {"1", "1.0", "m", "male"}:
        return 1
    if value in {"2", "2.0", "f", "female"}:
        return 2
    return None


def _normalized_chromosome(chrom):
    return str(chrom).strip().upper().removeprefix("CHR")



def check_variant_presence(alleles, gender, chrom):
    """Determine if genotype matches the inheritance mode, considering gender and chromosome."""
    if alleles is None:
        return None
    variant_count = sum(1 for allele in alleles if str(allele) == '1')
    chromosome = _normalized_chromosome(chrom)
    sex = _normalized_sex(gender)
    is_x = chromosome == 'X'
    is_y = chromosome == 'Y'

    if is_y:
        if sex == 1:
            return variant_count >= 1
        if sex == 2:
            return False
        return None

    if is_x:
        if sex == 1:
            return "Hemizygous" if variant_count >= 1 else False
        if sex == 2:
            if variant_count > 1:
                return "Homozygous"
            elif variant_count == 1:
                return "Heterozygous"
            else:
                return False
        return None

    # Autosomal (AD or AR)
    if variant_count > 1:
        return "Homozygous"
    elif variant_count == 1:
        return "Heterozygous"
    else:
        return False



def find_cosegregating_variants(vcf_file, ped_file, mode="both"):
    """Count relatives whose genotypes cosegregate with each variant.
    
    Args:
        vcf_file (str): Path to multi-family VCF file.
        ped_file (str): Path to pedigree file.
        mode (str): Compatibility selector accepted by the public CLI.
    Returns:
        dict: Per-variant segregation counts, including explicit male and female
            counts for sex-linked criteria.
    """

    # Load pedigree
    ped_df = pd.read_csv(ped_file, sep='\t', header=None, 
                         names=['FamilyID', 'IndividualID', 'PaternalID', 'MaternalID', 'Sex', 'Phenotype'])
    if mode not in {"dominant", "recessive", "both"}:
        raise ValueError(f"Unsupported cosegregation mode: {mode}")
    ped_df.loc[:, 'Sex'] = pd.to_numeric(ped_df['Sex'], errors='coerce').astype('Int64')
    ped_df.loc[:, 'Phenotype'] = pd.to_numeric(
        ped_df['Phenotype'], errors='coerce'
    ).astype('Int64')
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
            affected_male_segregation_count = 0
            unaffected_male_segregation_count = 0
            affected_female_segregation_count = 0
            unaffected_female_segregation_count = 0
            # Initialize results
            results = { "Affected_segregated_inds": 0,
                        "Unaffected_segregated_inds": 0,
                        "Affected_homozygous_inds": 0,
                        "Unaffected_nonhomo_inds": 0,
                        "Male_affected_segregated_inds": 0,
                        "Male_unaffected_segregated_inds": 0,
                        "Female_affected_segregated_inds": 0,
                        "Female_unaffected_segregated_inds": 0 }
    
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
                female_affected_gts = []
                for ind, gender in affected:
                    sex = _normalized_sex(gender)
                    if _normalized_chromosome(chrom) in {"X", "Y"} and sex is None:
                        continue
                    gt = record.samples[ind]['GT']
                    alleles = parse_gt(gt)
                    presence = check_variant_presence(alleles, gender, chrom)
                    if not presence:
                        all_families_pass = False
                        break
                    else:
                        affected_ind_gts.append(presence)
                        if sex == 1:
                            male_affected_gts.append(presence)
                        elif sex == 2:
                            female_affected_gts.append(presence)

                if not all_families_pass:
                    # Early exit if affected individuals do not carry the variant
                    break

                # Check unaffected individuals
                non_parental_controls_gts = []
                parental_controls_gts = []
                male_nonparental_control_gts = []
                female_nonparental_control_gts = []

                for ind, gender in unaffected:
                    sex = _normalized_sex(gender)
                    if _normalized_chromosome(chrom) in {"X", "Y"} and sex is None:
                        continue
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
                            if sex == 1:
                                male_nonparental_control_gts.append(presence)
                            elif sex == 2:
                                female_nonparental_control_gts.append(presence)

                if not all_families_pass:
                    # Early exit if any unaffected individual carries the variant in disallowed state
                    break

                control_gts = parental_controls_gts + non_parental_controls_gts

                if not any(control_gts):
                    affected_segregation_count += len(affected_ind_gts)
                    unaffected_segregation_count += len(non_parental_controls_gts)
                    affected_male_segregation_count += len(male_affected_gts)
                    unaffected_male_segregation_count += len(male_nonparental_control_gts)
                    affected_female_segregation_count += len(female_affected_gts)
                    unaffected_female_segregation_count += len(
                        female_nonparental_control_gts
                    )
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
                results["Male_unaffected_segregated_inds"] += unaffected_male_segregation_count
                results["Female_affected_segregated_inds"] += affected_female_segregation_count
                results["Female_unaffected_segregated_inds"] += unaffected_female_segregation_count
                per_var_stats[variant] = results
               
    return per_var_stats



if __name__ == "__main__":
    parser = ap.ArgumentParser()
    parser.add_argument("--vcf", type=str, required=True)
    parser.add_argument("--ped", type=str, required=True)
    parser.add_argument("--mode", type=str, required=True, choices=['dominant', 'recessive', 'both'])
    args = parser.parse_args()
    results = find_cosegregating_variants(args.vcf, args.ped, args.mode)
