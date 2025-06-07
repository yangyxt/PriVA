import pandas as pd
import numpy as np
from collections import defaultdict
from typing import Callable, Optional, List, Dict, Any, NamedTuple, Tuple
import multiprocessing
import os # For os.cpu_count()

logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)
console_handler=logging.StreamHandler()
console_handler.setLevel(logging.INFO)
formatter = logging.Formatter("%(levelname)s:%(asctime)s:%(funcName)s:%(lineno)s:%(message)s")
console_handler.setFormatter(formatter)
logger.addHandler(console_handler)

# --- START: determine_cis_trans and its helpers ---
# This section is identical to the previous response.
# It includes ParsedGenotype, parse_genotype_string,
# get_cis_trans_from_haplotypes, can_parent_transmit_haplotype,
# and determine_cis_trans (with male X variants as "in-cis").
class ParsedGenotype(NamedTuple):
    allele1: Optional[str]
    allele2: Optional[str]
    is_phased: bool
    is_fully_known: bool
    has_alt_allele: bool

def parse_genotype_string(gt_str: Optional[str], gender: str, chromosome: str) -> ParsedGenotype:
    if gt_str is None:
        return ParsedGenotype(None, None, False, False, False)

    chrom_upper = str(chromosome).upper() 
    is_male_x = gender.lower() == "male" and (chrom_upper == "X" or chrom_upper == "CHRX")
    
    phased_separator = '|'
    unphased_separator = '/'
    is_phased = phased_separator in gt_str
    a1_val: Optional[str] = None
    a2_val: Optional[str] = None

    if gt_str == ".": 
        if is_male_x: 
            a1_val, a2_val = '0', '0' 
            return ParsedGenotype(a1_val, a2_val, True, True, a1_val == '1')
        else: 
            return ParsedGenotype(None, None, False, False, False) 
    
    if gt_str == "./.":
         return ParsedGenotype('0', '0', False, True, False)

    if is_male_x:
        allele_val_str = None
        if phased_separator in gt_str or unphased_separator in gt_str:
            parts = gt_str.split(phased_separator if is_phased else unphased_separator)
            if len(parts) == 2 and parts[0] == parts[1] and parts[0] in ['0', '1', '.']:
                allele_val_str = parts[0]
            elif len(parts) == 1 and parts[0] in ['0', '1', '.']:
                 allele_val_str = parts[0]
            else: 
                return ParsedGenotype(None, None, False, False, False)
        elif gt_str in ['0', '1', '.']:
            allele_val_str = gt_str
        else:
            return ParsedGenotype(None, None, False, False, False)

        if allele_val_str == '.':
             a1_val, a2_val = '0', '0' 
        else:
            a1_val, a2_val = allele_val_str, allele_val_str
        
        is_fully_known = True 
        has_alt = (a1_val == '1') 
        return ParsedGenotype(a1_val, a2_val, True, is_fully_known, has_alt)

    if not (is_phased or unphased_separator in gt_str):
        return ParsedGenotype(None, None, False, False, False)

    parts = gt_str.split(phased_separator if is_phased else unphased_separator)
    if len(parts) != 2:
        return ParsedGenotype(None, None, False, False, False)
    
    s_a1, s_a2 = parts[0], parts[1]

    if s_a1 == '.' and s_a2 == '.':
        a1_val, a2_val = '0', '0'
    elif s_a1 == '.':
        a1_val = '0'
        a2_val = s_a2 if s_a2 in ['0','1'] else 'invalid_allele'
    elif s_a2 == '.':
        a1_val = s_a1 if s_a1 in ['0','1'] else 'invalid_allele'
        a2_val = '0'
    else:
        a1_val = s_a1 if s_a1 in ['0','1'] else 'invalid_allele'
        a2_val = s_a2 if s_a2 in ['0','1'] else 'invalid_allele'

    if a1_val == 'invalid_allele' or a2_val == 'invalid_allele':
        return ParsedGenotype(None, None, is_phased, False, False)

    is_fully_known = True
    has_alt = (a1_val == '1') or (a2_val == '1')
    return ParsedGenotype(a1_val, a2_val, is_phased, is_fully_known, has_alt)

def get_cis_trans_from_haplotypes(hap1: Tuple[str, str], 
                                  hap2: Tuple[str, str]) -> str:
    h1_v1, h1_v2 = hap1
    h2_v1, h2_v2 = hap2
    cis_on_hap1 = (h1_v1 == '1' and h1_v2 == '1')
    cis_on_hap2 = (h2_v1 == '1' and h2_v2 == '1')
    if cis_on_hap1 or cis_on_hap2: return "in-cis"
    trans_config1 = (h1_v1 == '1' and h2_v2 == '1') 
    trans_config2 = (h2_v1 == '1' and h1_v2 == '1')
    if trans_config1 or trans_config2: return "in-trans"
    return str(np.nan)

def can_parent_transmit_haplotype(
    parent_parsed_v1: Optional[ParsedGenotype], parent_parsed_v2: Optional[ParsedGenotype],
    target_haplotype: Tuple[str, str]) -> bool:
    v1_parent_gt_unknown = parent_parsed_v1 is None or \
                           (parent_parsed_v1.allele1 is None and parent_parsed_v1.allele2 is None)
    v2_parent_gt_unknown = parent_parsed_v2 is None or \
                           (parent_parsed_v2.allele1 is None and parent_parsed_v2.allele2 is None)
    if v1_parent_gt_unknown and v2_parent_gt_unknown: return True
    p_v1_a1 = parent_parsed_v1.allele1 if parent_parsed_v1 else None
    p_v1_a2 = parent_parsed_v1.allele2 if parent_parsed_v1 else None
    p_v2_a1 = parent_parsed_v2.allele1 if parent_parsed_v2 else None
    p_v2_a2 = parent_parsed_v2.allele2 if parent_parsed_v2 else None
    target_v1_allele, target_v2_allele = target_haplotype
    def check_parent_hap_compatibility(p_h_a1, p_h_a2):
        v1_match = (p_h_a1 is None or p_h_a1 == target_v1_allele)
        v2_match = (p_h_a2 is None or p_h_a2 == target_v2_allele)
        return v1_match and v2_match
    if parent_parsed_v1 and parent_parsed_v1.is_phased and \
       parent_parsed_v2 and parent_parsed_v2.is_phased:
        if check_parent_hap_compatibility(p_v1_a1, p_v2_a1): return True
        if check_parent_hap_compatibility(p_v1_a2, p_v2_a2): return True
        return False
    if p_v1_a1 is not None and p_v1_a1 == p_v1_a2: 
        if p_v1_a1 != target_v1_allele: return False
        if check_parent_hap_compatibility(p_v1_a1, p_v2_a1): return True
        if check_parent_hap_compatibility(p_v1_a1, p_v2_a2): return True
        return False
    if p_v2_a1 is not None and p_v2_a1 == p_v2_a2: 
        if p_v2_a1 != target_v2_allele: return False
        if check_parent_hap_compatibility(p_v1_a1, p_v2_a1): return True
        if check_parent_hap_compatibility(p_v1_a2, p_v2_a1): return True
        return False
    v1_can_be_supplied = (p_v1_a1 is None and p_v1_a2 is None) or \
                         (p_v1_a1 == target_v1_allele) or (p_v1_a2 == target_v1_allele)
    v2_can_be_supplied = (p_v2_a1 is None and p_v2_a2 is None) or \
                         (p_v2_a1 == target_v2_allele) or (p_v2_a2 == target_v2_allele)
    return v1_can_be_supplied and v2_can_be_supplied

def determine_cis_trans(
    sample_var1_gt_str: Optional[str], sample_var2_gt_str: Optional[str], 
    gender: str, chromosome: str,
    father_var1_gt_str: Optional[str] = None, father_var2_gt_str: Optional[str] = None,
    mother_var1_gt_str: Optional[str] = None, mother_var2_gt_str: Optional[str] = None) -> str:
    str_nan = str(np.nan) 
    chrom_upper = str(chromosome).upper()
    sample_p_v1 = parse_genotype_string(sample_var1_gt_str, gender, chromosome)
    sample_p_v2 = parse_genotype_string(sample_var2_gt_str, gender, chromosome)
    if sample_p_v1.allele1 is None or sample_p_v1.allele2 is None or \
       sample_p_v2.allele1 is None or sample_p_v2.allele2 is None: return str_nan
    if not sample_p_v1.has_alt_allele or not sample_p_v2.has_alt_allele: return str_nan 
    if gender.lower() == "male" and (chrom_upper == "X" or chrom_upper == "CHRX"): return "in-cis" 
    s_v1_a1, s_v1_a2 = str(sample_p_v1.allele1), str(sample_p_v1.allele2)
    s_v2_a1, s_v2_a2 = str(sample_p_v2.allele1), str(sample_p_v2.allele2)
    if sample_p_v1.is_phased and sample_p_v2.is_phased:
        return get_cis_trans_from_haplotypes((s_v1_a1, s_v2_a1), (s_v1_a2, s_v2_a2))
    father_p_v1 = parse_genotype_string(father_var1_gt_str, "male", chromosome)
    father_p_v2 = parse_genotype_string(father_var2_gt_str, "male", chromosome)
    mother_p_v1 = parse_genotype_string(mother_var1_gt_str, "female", chromosome)
    mother_p_v2 = parse_genotype_string(mother_var2_gt_str, "female", chromosome)
    child_hap1_config1, child_hap2_config1 = (s_v1_a1, s_v2_a1), (s_v1_a2, s_v2_a2)
    child_hap1_config2, child_hap2_config2 = (s_v1_a1, s_v2_a2), (s_v1_a2, s_v2_a1)
    consistency_1a = can_parent_transmit_haplotype(father_p_v1, father_p_v2, child_hap1_config1) and \
                     can_parent_transmit_haplotype(mother_p_v1, mother_p_v2, child_hap2_config1)
    consistency_1b = can_parent_transmit_haplotype(father_p_v1, father_p_v2, child_hap2_config1) and \
                     can_parent_transmit_haplotype(mother_p_v1, mother_p_v2, child_hap1_config1)
    config1_is_possible = consistency_1a or consistency_1b
    consistency_2a = can_parent_transmit_haplotype(father_p_v1, father_p_v2, child_hap1_config2) and \
                     can_parent_transmit_haplotype(mother_p_v1, mother_p_v2, child_hap2_config2)
    consistency_2b = can_parent_transmit_haplotype(father_p_v1, father_p_v2, child_hap2_config2) and \
                     can_parent_transmit_haplotype(mother_p_v1, mother_p_v2, child_hap1_config2)
    config2_is_possible = consistency_2a or consistency_2b
    status_config1_obj, status_config2_obj = np.nan, np.nan
    if config1_is_possible: status_config1_obj = get_cis_trans_from_haplotypes(child_hap1_config1, child_hap2_config1)
    if config2_is_possible: status_config2_obj = get_cis_trans_from_haplotypes(child_hap1_config2, child_hap2_config2)
    status_config1_str, status_config2_str = str(status_config1_obj), str(status_config2_obj)
    if config1_is_possible and not config2_is_possible: return status_config1_str
    if not config1_is_possible and config2_is_possible: return status_config2_str
    if config1_is_possible and config2_is_possible:
        return status_config1_str if status_config1_obj is not np.nan and status_config1_obj == status_config2_obj else str_nan
    return str_nan
# --- END: determine_cis_trans and its helpers ---

# --- Worker function for multiprocessing ---
# This function must be defined at the top level of a module to be picklable by multiprocessing
def _process_gene_for_patient_worker(
    # gene_name: str, # Not strictly needed by worker if just returning sets of g_var_keys
    variants_in_gene_list: List[Dict[str, Any]], 
    patient_sex: str,
    str_nan_val: str
    ) -> Tuple[set, set]:
    """
    Processes variants within a single gene for a single patient to find
    genomic variant keys that are in cis or in trans with a pathogenic variant.
    """
    gene_g_vars_cis_with_patho = set()
    gene_g_vars_trans_with_patho = set()

    if len(variants_in_gene_list) < 2:
        return gene_g_vars_cis_with_patho, gene_g_vars_trans_with_patho

    for i in range(len(variants_in_gene_list)):
        for j in range(i + 1, len(variants_in_gene_list)):
            var_A_info = variants_in_gene_list[i]
            var_B_info = variants_in_gene_list[j]

            path_A = var_A_info["is_pathogenic"]
            path_B = var_B_info["is_pathogenic"]

            if not (path_A or path_B): # At least one must be pathogenic
                continue
            
            chromosome = var_A_info["chrom"] 
            
            cis_trans_status = determine_cis_trans( # This function is globally defined
                sample_var1_gt_str=var_A_info["proband_gt"],
                sample_var2_gt_str=var_B_info["proband_gt"],
                gender=patient_sex,
                chromosome=chromosome,
                father_var1_gt_str=var_A_info["father_gt"],
                father_var2_gt_str=var_B_info["father_gt"], # Var B's father GT for var B
                mother_var1_gt_str=var_A_info["mother_gt"],
                mother_var2_gt_str=var_B_info["mother_gt"]  # Var B's mother GT for var B
            )

            if cis_trans_status != str_nan_val:
                g_var_key_A = var_A_info["id"]
                g_var_key_B = var_B_info["id"]
                if cis_trans_status == "in-cis":
                    if path_A: gene_g_vars_cis_with_patho.add(g_var_key_B)
                    if path_B: gene_g_vars_cis_with_patho.add(g_var_key_A)
                elif cis_trans_status == "in-trans":
                    if path_A: gene_g_vars_trans_with_patho.add(g_var_key_B)
                    if path_B: gene_g_vars_trans_with_patho.add(g_var_key_A)
    
    return gene_g_vars_cis_with_patho, gene_g_vars_trans_with_patho

# --- Main Batch Annotation Function ---
def batch_annotate_cis_trans_from_table(
    variants_df: pd.DataFrame,
    is_pathogenic_array: List[bool], 
    pedigree_df: pd.DataFrame,
    num_processes: Optional[int] = None, # For multiprocessing
    # --- Column name parameters for variants_df ---
    chrom_col: str = "Chromosome",
    pos_col: str = "Position",
    ref_col: str = "Reference_Allele",
    alt_col: str = "Alternate_Allele",
    gene_col: str = "Gene_Symbol",
    # --- Column name parameters for pedigree_df ---
    ped_sample_id_col: str = "IndividualID",
    ped_paternal_id_col: str = "PaternalID",
    ped_maternal_id_col: str = "MaternalID",
    ped_sex_col: str = "Sex",
    ped_phenotype_col: str = "Phenotype", 
    ped_patient_value: Any = 2) -> Tuple[np.ndarray, np.ndarray]:
    """
    Processes a variant table to find pairs of variants within the same gene for "patients",
    determines their phase using multiprocessing for per-gene tasks per patient, 
    and returns two integer arrays counting for how many patients each variant-transcript 
    row is in cis or in trans with another pathogenic variant.
    """
    str_nan = str(np.nan)

    # --- 1. Data Preparation ---
    if len(variants_df) != len(is_pathogenic_array):
        raise ValueError("Length of variants_df and is_pathogenic_array must be the same.")
    
    df_proc = variants_df.copy()
    df_proc["_is_pathogenic_transcript"] = is_pathogenic_array
    df_proc["_original_row_index"] = np.arange(len(df_proc))
    
    df_proc[chrom_col] = df_proc[chrom_col].astype(str)
    df_proc[pos_col] = df_proc[pos_col].astype(str) 
    df_proc[ref_col] = df_proc[ref_col].astype(str)
    df_proc[alt_col] = df_proc[alt_col].astype(str)
    df_proc[gene_col] = df_proc[gene_col].astype(str)

    df_proc["_genomic_variant_key"] = df_proc[chrom_col] + "_" + \
                                     df_proc[pos_col] + "_" + \
                                     df_proc[ref_col] + "_" + \
                                     df_proc[alt_col]

    genomic_variant_pathogenicity = df_proc.groupby("_genomic_variant_key")["_is_pathogenic_transcript"].any().rename("_is_g_var_pathogenic")
    
    unique_genomic_variants_grouped = df_proc.groupby("_genomic_variant_key").agg(
        **{ 
            chrom_col: pd.NamedAgg(column=chrom_col, aggfunc='first'),
            pos_col: pd.NamedAgg(column=pos_col, aggfunc='first'),
            ref_col: pd.NamedAgg(column=ref_col, aggfunc='first'),
            alt_col: pd.NamedAgg(column=alt_col, aggfunc='first'),
            gene_col: pd.NamedAgg(column=gene_col, aggfunc=lambda x: set(g for g in x.dropna().unique() if g and pd.notna(g)))
        }
    ).reset_index()
    unique_genomic_variants = unique_genomic_variants_grouped.merge(genomic_variant_pathogenicity, on="_genomic_variant_key", how="left")
    unique_genomic_variants[pos_col] = unique_genomic_variants[pos_col].astype(int)

    g_var_key_to_original_indices = df_proc.groupby("_genomic_variant_key")['_original_row_index'].apply(lambda x: x.to_numpy(copy=False)).to_dict()

    patients_in_cis_count_arr = np.zeros(len(df_proc), dtype=int)
    patients_in_trans_count_arr = np.zeros(len(df_proc), dtype=int)

    if ped_phenotype_col not in pedigree_df.columns:
        raise ValueError(f"Phenotype column '{ped_phenotype_col}' not found in pedigree DataFrame.")
    
    patient_ped_rows_df = pedigree_df[pedigree_df[ped_phenotype_col].astype(int) == ped_patient_value]
    if patient_ped_rows_df.empty:
        logger.warning("Warning: No individuals identified as 'patients'. Returning zero-filled count arrays.")
        return patients_in_cis_count_arr, patients_in_trans_count_arr

    actual_num_processes = 1 if num_processes is None else num_processes
    logger.info(f"Using {actual_num_processes} processes for per-gene analysis per patient.")


    for _, patient_ped_row in patient_ped_rows_df.iterrows():
        proband_id = str(patient_ped_row[ped_sample_id_col])
        # ... (proband_id, father_id, mother_id, proband_sex extraction - same as before) ...
        father_id_val = patient_ped_row.get(ped_paternal_id_col)
        mother_id_val = patient_ped_row.get(ped_maternal_id_col)
        father_id = str(father_id_val) if pd.notna(father_id_val) and str(father_id_val) != "0" and father_id_val != 0 else None
        mother_id = str(mother_id_val) if pd.notna(mother_id_val) and str(mother_id_val) != "0" and mother_id_val != 0 else None
        proband_sex_val = patient_ped_row.get(ped_sex_col)
        proband_sex = ""
        if isinstance(proband_sex_val, (int, float)):
            if proband_sex_val == 1: proband_sex = "male"
            elif proband_sex_val == 2: proband_sex = "female"
            else: proband_sex = "unknown"
        elif isinstance(proband_sex_val, str):
            proband_sex_val_lower = proband_sex_val.lower()
            if proband_sex_val_lower in ["male", "female"]: proband_sex = proband_sex_val_lower
            else: proband_sex = "unknown"
        else: proband_sex = "unknown"

        if proband_sex == "unknown":
            logger.warning(f"Warning: Patient {proband_id} has unknown sex. Skipping.")
            continue
        if proband_id not in df_proc.columns:
            logger.debug(f"Warning: GT column for patient {proband_id} not found. Skipping.")
            continue
        if father_id and father_id not in df_proc.columns: father_id = None
        if mother_id and mother_id not in df_proc.columns: mother_id = None

        logger.info(f"Processing patient {proband_id} with sex {proband_sex} and father {father_id} and mother {mother_id}")

        proband_variants_by_gene = defaultdict(list)
        for _, g_var_row in unique_genomic_variants.iterrows():
            g_var_key = g_var_row["_genomic_variant_key"]
            first_original_row = df_proc[df_proc["_genomic_variant_key"] == g_var_key].iloc[0]
            variant_info = {
                "id": g_var_key, "chrom": str(g_var_row[chrom_col]), "pos": int(g_var_row[pos_col]),
                "ref": str(g_var_row[ref_col]), "alt": str(g_var_row[alt_col]),
                "is_pathogenic": bool(g_var_row["_is_g_var_pathogenic"]),
                "proband_gt": str(first_original_row.get(proband_id, "./.")),
                "father_gt": str(first_original_row.get(father_id, "./.")) if father_id else None,
                "mother_gt": str(first_original_row.get(mother_id, "./.")) if mother_id else None
            }
            for gene_name in g_var_row[gene_col]:
                if gene_name and pd.notna(gene_name):
                    proband_variants_by_gene[str(gene_name)].append(variant_info)
        
        # Prepare tasks for multiprocessing pool for this patient
        pool_tasks = []
        for gene_name, variants_in_gene_list in proband_variants_by_gene.items():
            if len(variants_in_gene_list) >= 2:
                # Pass variants_in_gene_list, patient_sex, str_nan
                pool_tasks.append((variants_in_gene_list, proband_sex, str_nan))
        
        g_vars_cis_with_patho_this_patient = set()
        g_vars_trans_with_patho_this_patient = set()

        if pool_tasks:
            if actual_num_processes > 1 and len(pool_tasks) > 1 : # Only use pool if multiple tasks and processes
                with multiprocessing.Pool(processes=actual_num_processes) as pool:
                    results_from_pool = pool.starmap(_process_gene_for_patient_worker, pool_tasks)
                
                for gene_cis_set, gene_trans_set in results_from_pool:
                    g_vars_cis_with_patho_this_patient.update(gene_cis_set)
                    g_vars_trans_with_patho_this_patient.update(gene_trans_set)
            else: # Process serially if only 1 process or 1 task
                for task_args in pool_tasks:
                    gene_cis_set, gene_trans_set = _process_gene_for_patient_worker(*task_args)
                    g_vars_cis_with_patho_this_patient.update(gene_cis_set)
                    g_vars_trans_with_patho_this_patient.update(gene_trans_set)

        # Update global count arrays based on this patient's findings
        for g_key in g_vars_cis_with_patho_this_patient:
            indices = g_var_key_to_original_indices.get(g_key)
            if indices is not None: patients_in_cis_count_arr[indices] += 1
        
        for g_key in g_vars_trans_with_patho_this_patient:
            indices = g_var_key_to_original_indices.get(g_key)
            if indices is not None: patients_in_trans_count_arr[indices] += 1
            
    return patients_in_cis_count_arr, patients_in_trans_count_arr