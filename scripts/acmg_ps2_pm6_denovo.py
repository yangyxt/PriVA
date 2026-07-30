#!/usr/bin/env python3
"""PS2 and PM6 -- the variant arose new in the proband.

    PS2  confirmed de novo: present in the proband, absent from both parents,
         and parentage confirmed
    PM6  assumed de novo: only one parent available, and the allele is absent
         from gnomAD

identify_fam_members resolves the proband, parents and siblings from the
pedigree table; determine_denovo_per_row compares their genotypes.
"""

import logging
import pandas as pd
import numpy as np
from typing import Tuple


logger = logging.getLogger(__name__)


def identify_fam_members(ped_df: pd.DataFrame, fam_name: str) -> pd.DataFrame:
    """
    Identify proband/parents/siblings for a given family from a PED table.

    IMPORTANT:
    - Father/mother must come from the *proband's* PaternalID/MaternalID fields.
    """
    fam_ped_df = ped_df.loc[ped_df['#FamilyID'] == fam_name, :].copy()
    if fam_ped_df.empty:
        raise ValueError(f"Family {fam_name} not found in pedigree table")

    # Normalize types (PEDs sometimes store these as strings)
    fam_ped_df.loc[:, "Phenotype"] = pd.to_numeric(fam_ped_df["Phenotype"], errors="coerce").astype("Int64")
    fam_ped_df.loc[:, "Sex"] = pd.to_numeric(fam_ped_df["Sex"], errors="coerce")

    patients = fam_ped_df.loc[fam_ped_df["Phenotype"] == 2, "IndividualID"].tolist()
    if not patients:
        raise ValueError(f"No affected individual (Phenotype==2) found for family {fam_name}")
    if len(patients) > 1:
        logger.warning(f"Family {fam_name} has multiple affected individuals (Phenotype==2): {patients}. Using the first as proband.")
    proband = patients[0]

    proband_row = fam_ped_df.loc[fam_ped_df["IndividualID"] == proband, :].iloc[0]

    na_values = {"0", 0, ".", "", "NA", "NaN", "nan", -9, "-9"}
    father = proband_row.get("PaternalID", None)
    mother = proband_row.get("MaternalID", None)
    father = None if (pd.isna(father) or father in na_values) else str(father)
    mother = None if (pd.isna(mother) or mother in na_values) else str(mother)

    father_pheno = None
    if father and (father in set(fam_ped_df["IndividualID"].astype(str))):
        father_pheno = fam_ped_df.loc[fam_ped_df["IndividualID"].astype(str) == father, "Phenotype"].iloc[0]
    mother_pheno = None
    if mother and (mother in set(fam_ped_df["IndividualID"].astype(str))):
        mother_pheno = fam_ped_df.loc[fam_ped_df["IndividualID"].astype(str) == mother, "Phenotype"].iloc[0]

    # Process siblings: all family members excluding proband and parents (if present)
    exist_mems = {m for m in [father, mother, proband] if m is not None}
    sib_df = fam_ped_df.loc[~fam_ped_df["IndividualID"].astype(str).isin(exist_mems), ["IndividualID", "Phenotype"]]
    sib_info = dict(zip(sib_df["IndividualID"].astype(str).tolist(), sib_df["Phenotype"].tolist()))

    return (proband, 2), (father, father_pheno), (mother, mother_pheno), sib_info


def determine_denovo_per_row(row: dict, ped_info: dict, ped_df: pd.DataFrame) -> bool:
    # Input row is a dictionary converted from a pandas Series
    # Determine whether the variant is a denovo mutation in the proband
    proband_info, father_info, mother_info, sib_info = ped_info
    proband, proband_pheno = proband_info
    father, father_pheno = father_info
    mother, mother_pheno = mother_info

    # Coerce NaN genotype values to '.' so that `"1" in gt` never hits a float
    for sample_id in (father, mother, proband):
        val = row.get(sample_id)
        if val is not None and not isinstance(val, str):
            row[sample_id] = '.'

    proband_sex = ped_df.loc[ped_df['IndividualID'] == proband, 'Sex'].values[0]
    # First determine when proband is Male.
    if proband_sex == "1" or proband_sex == 1:
        if row['chrom'] == "chrX":
            # Variant on chrX are hemizygous in males
            if "1" in row.get(mother, '.'):
                return False
            elif not "." in row.get(mother, '.'):
                return "PS2"
        elif row['chrom'] == "chrY":
            # Variant on chrY are hemizygous in males
            if "1" in row.get(father, '.'):
                return False
            elif not "." in row.get(father, '.'):
                return "PS2"
        elif row['chrom'] == "chrM":
            # Variant on chrM are hemizygous in males
            if "1" in row.get(mother, '.'):
                return False
            elif not "." in row.get(mother, '.'):
                return "PS2"
        else:
            if "1" in row.get(father, '.') or "1" in row.get(mother, '.'):
                return False
            elif (not "." in row.get(father, '.')) and (not "." in row.get(mother, '.')):
                return "PS2"
            elif row.get('gnomAD_joint_AF') in [0, np.nan] and ((not "." in row.get(father, '.')) or (not "." in row.get(mother, '.'))):
                return "PM6"
            else:
                return False
    else:
        if row['chrom'] == "chrX":
            if "1" in row.get(father, '.') or "1" in row.get(mother, '.'):
                return False
            elif (not "." in row.get(father, '.')) and (not "." in row.get(mother, '.')):
                return "PS2"
            elif row.get('gnomAD_joint_AF') in [0, np.nan] and ((not "." in row.get(father, '.')) or (not "." in row.get(mother, '.'))):
                return "PM6"
            else:
                return False
        elif row['chrom'] == "chrY":
            if "1" in row.get(father, '.'):
                return False
            elif (not "." in row.get(father, '.')):
                return "PS2"
        elif row['chrom'] == "chrM":
            if "1" in row.get(mother, '.'):
                return False
            elif (not "." in row.get(mother, '.')):
                return "PS2"
        else:
            if "1" in row.get(father, '.') or "1" in row.get(mother, '.'):
                return False
            elif (not "." in row.get(father, '.')) and (not "." in row.get(mother, '.')):
                return "PS2"
            elif row.get('gnomAD_joint_AF') in [0, np.nan] and ((not "." in row.get(father, '.')) or (not "." in row.get(mother, '.'))):
                return "PM6"
            else:
                return False


def PS2_PM6_criteria(df: pd.DataFrame,
                     ped_df: pd.DataFrame,
                     fam_name: str,
                     threads: int = 10) -> Tuple[np.ndarray, np.ndarray]:
    """
    PS2: confirmed denovo mutation in the proband
    PM6: Assumed denovo mutation in the proband (if only have one parent and the PAF is absent from gnomAD)
    Currently, we cant implement the variable strength level because the phenotype-gene relationship cant be evaluated by purely bioinformatic means
    """
    logger.info(f"Running determine_denovo using pandas apply method")
    proband_info, father_info, mother_info, sib_info = identify_fam_members(ped_df, fam_name)
    proband, _ = proband_info
    father, _ = father_info
    mother, _ = mother_info

    # Get subset of pedigree for this family
    ped_subset = ped_df.loc[ped_df['#FamilyID'] == fam_name, :]

    # Select only necessary columns to reduce memory usage
    # These columns are needed for the denovo determination
    necessary_cols = ['chrom', 'pos', 'ref', 'alt', 'gnomAD_joint_AF']

    # Add columns for family members if they exist
    if proband in df.columns:
        necessary_cols.append(proband)
    if father in df.columns:
        necessary_cols.append(father)
    if mother in df.columns:
        necessary_cols.append(mother)

    # Create a unique key for each variant
    df['variant_key'] = df['chrom'] + '_' + df['pos'].astype(str) + '_' + df['ref'] + '_' + df['alt']

    # Keep track of original index
    original_index = df.index.copy()

    # Select necessary columns and drop duplicates
    cols_to_use = [col for col in necessary_cols if col in df.columns]
    unique_df = df[cols_to_use + ['variant_key']].drop_duplicates(subset=['variant_key'])
    unique_df['gnomAD_joint_AF'] = unique_df['gnomAD_joint_AF'].fillna(0)

    logger.info(f"Processing {len(unique_df)} unique variants (reduced from {len(df)} total variants)")

    # Create ped_info tuple
    ped_info = (proband_info, father_info, mother_info, sib_info)

    # Apply the function to each unique row
    unique_results = unique_df.apply(
        lambda row: determine_denovo_per_row(row.to_dict(), ped_info, ped_subset),
        axis=1
    )

    # Create a mapping from variant keys to results
    result_map = dict(zip(unique_df['variant_key'], unique_results))

    # Map results back to original dataframe
    results = np.array([result_map.get(key, False) for key in df['variant_key']])

    # Clean up temporary column
    df.drop('variant_key', axis=1, inplace=True)

    # Convert results to final form
    ps2_criteria = results == "PS2"
    pm6_criteria = results == "PM6"

    logger.info(f"Found {ps2_criteria.sum()} variants with PS2 criteria (confirmed de novo)")
    logger.info(f"Found {pm6_criteria.sum()} variants with PM6 criteria (assumed de novo)")

    ps2_array = np.zeros(len(df), dtype=int)
    pm6_array = np.zeros(len(df), dtype=int)
    ps2_array[ps2_criteria] = 3
    pm6_array[pm6_criteria] = 2

    return ps2_array, pm6_array
